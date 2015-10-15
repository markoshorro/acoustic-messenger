function [Xhat, psd, const, eyed] = receiver(tout,fc)
	%% RECEIVER FUNCTION
    % Group 13
    % Introduction to Communication Engineering. September 2015 
    %
    % Implementation of the receiver.
    % INPUT:  tout  - type double and defines the time during while
    %                 waiting
    %         fc    - carrier frequency
    %
    % OUTPUT: Xhat  - Xhat is a row-vector of 432 bits (0, 1) of the
    %                 received information.
    %         psd   - psd is a structure with two fields: psd.p is a row
    %                 vector of the PSD of the received signal calculated
    %                 after down-conversion in dB. The vector psd.p is normalized
    %                 that normalized so that maximum value is 0 dB. psd.f is a
    %                 row-vector of corresponding frequencies in Hz
    %         const - row-vector of complex samples after downsampling that 
    %                 correspond to information par
    %         eyed  - eyed.fsfd is an integer and stands for a number of 
    %                 samples per symbol; eyed.r is a complex row-vector of
    %                 information part of the received signal after matched
    %                 filtering and timing synchronization (properly cut so
    %                 that MATLAB eyediagram(eyed.r, eyed.fsfd) plots a 
    %                 desired eye diagram).
    %
	%% Some parameters
    run('../parameters.m');
    
    %% Recorder declarations
    channels = 1;       % Mono
    recordBits = 16;    % More precision
    % Creating recording object
    recording = audiorecorder(fs, recordBits, channels);   
    
    %% Recording Audio
    % Create the barker pulse train
    [si,~] = rtrcpuls(rollOff, Tau, fs, span);
    symbolsBarker = constBPSK(symbBarker);
    pulseBarker = conv(upsample(symbolsBarker, sps), si);
    
    %% Declarations
    % Check if Barker was found
    foundBarker = 0;
    % Difference between peaks on a window correlated
    thresholdDiff = 20;
    % Threshold 
    thresholdPeak = 20;
    % Pause between iterations on the loop
    timePreamble = .2;
    % Time to record all the data after finding the Barker sequence
    timeData = 1;
    % Index for the window
    idxWin = 1;
    % Start Recording
    record(recording);
    tic;
    
    % Recording loop
    while (true)
        % Record 1 second (entire packet is 0.99s)
        pause(timePreamble);
            
        % Fetch Data
        inputSound = getaudiodata(recording, 'double');
        tmpWindow = inputSound(idxWin:end);
        
        % Converting window to baseband
        t = ((1:length(tmpWindow))/fs).';
        inputData = tmpWindow.*(exp(-1i*2*pi*fc*t));
        
        % Convolution with the pulse
        yt = conv(si, inputData);
        
        % If Barker sequence found, then we do not need to keep looking for
        % more data since we waited on the previous iteration, then we skip
        % the loop and continue
        if (foundBarker)
            break
        end
        
        % Updating the window index
        prevIdxWin = idxWin;
        idxWin = length(inputSound);

        % Cross correlation with the barker pulse to find the sequence
        corr = conv(yt, fliplr(pulseBarker));
        
        % Check for correlation with barker code
        [peakV, idxPeak] = max(abs(corr));
        corr(idxPeak-(floor(nBarker/2))*sps:...
            idxPeak+(floor(nBarker/2))*sps) = [];
        
        % Check for difference between barker code and packet
        [peakV2, ~] = max(abs(corr));
        diff = peakV - peakV2;

        % If there is a peak, then we found the Barker sequence
        if (peakV > thresholdPeak) && (diff > thresholdDiff)
            idxPreamble = idxPeak;
            foundBarker = 1;
            pause(timeData);
            idxWin = prevIdxWin;
        end
            
        % If the preamble has not been found in tout sec, return
        if (tout < toc)
            disp('I found nothing')
            Xhat=[]; psd=[]; const=[]; eyed=[];            
            return
        end
    end
    stop(recording);

    % Calculate the delay so that the start of the packet is known    
    delayHat = idxPreamble - length(pulseBarker) + sps*span + 1;
    
    % Remove all bits before the matching in convolution
    uncutPacket = yt(delayHat:end);
    uncutPacket = uncutPacket(1:(Ns+nBarker+nGuard)*sps);
    
    %% Decision: correct for QPSK
    % Packet downsampled
    downPacket = downsample(uncutPacket, sps);
    % Barker samples
    constBarker = downPacket(1:nBarker);
    % Data samples
    constPack = downPacket(nBarker+nGuard+1:nBarker+nGuard+Ns);
    
    % Calculating Phaseshift: Barker way   
    tmpOnes = constBarker(barker == 1);
    tmpOnes = [mean(real(tmpOnes)) mean(imag(tmpOnes))];
    tmpOnes = tmpOnes(:,1) + tmpOnes(:,2)*1i;
    phaseShift = wrapTo2Pi(angle(tmpOnes))-5*pi/4;
    
    %% Const output
    const = constPack*exp(-1i*phaseShift);
    
    %% Matching samples
    symbolsRec = const;
    symbolsRec(real(symbolsRec) > 0 & imag(symbolsRec) > 0) = 0;
    symbolsRec(real(symbolsRec) > 0 & imag(symbolsRec) < 0) = 1;
    symbolsRec(real(symbolsRec) < 0 & imag(symbolsRec) > 0) = 2;
    symbolsRec(real(symbolsRec) < 0 & imag(symbolsRec) < 0) = 3;
    
    %% Bits matched
    bitsGroup = de2bi(symbolsRec, m, 'left-msb');
        
    %% Outputs
    Xhat = reshape(bitsGroup.',[1,m*length(bitsGroup)]);

    % Normalization process
    tempConst = zeros(1,length(const));
    maxNorm = 0; 
    
    % Looking for the maximun normalized value
    for j=1:length(const)
        if(norm(const(j))>maxNorm)
            maxNorm = norm(const(j));
        end
    end

    % Dividing all values with the maximum value found before
    for i=1:length(tempConst)    
     tempConst(i)=const(i)./maxNorm;
    end
    
    % Output constellation normalized
	const = tempConst;    
      
    %% PSD output
    % Taking the passband signal
    psdData = tmpWindow;
    % Selecting passband singal related with the data
    psdData = psdData((nBarker+nGuard)*sps:(Ns+nBarker+nGuard)*sps);
    % Getting the PSD via Welch's method
    [pwData,pwVec] = pwelch(psdData,[],[],[],fs,'twosided');
    % Length of each segment in pwVec
    stepPw = (length(pwVec)/fs); 
    % Conversion to dB and normalization (to 0 db)
    pwData = 10*log10(pwData./max(pwData((round(stepPw*fc)-...
        round(stepPw*500):round(stepPw*fc)+round(stepPw*500)))));
    % Output PSD
    psd = struct('p',pwData,'f',pwVec-fc);
    
    %% Eyed output
    correctedPacket = uncutPacket((1+nBarker+nGuard)*sps:...
        sps*(nBarker+nGuard+Ns-1))*exp(-1i*phaseShift);
    eyed = struct('r',correctedPacket,'fsfd',sps);

end
