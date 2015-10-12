%function [Xhat, psd, const, eyed] = receiver(tout,fc)
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
    run('../parameters.m')
%    fc = 5000;
%    tout=inf;
    
    %% Recorder declarations
    channels = 1;       % Mono
    recordBits = 16;    % More precision
    
    % Creating recording object
    recording = audiorecorder(fs, recordBits, channels);   
    
    %% Low pass filter
    
    order = 2^10;
    lpf = fir1(order, Wn);
    
    %% Recording Audio

    % Create the barker pulse train
    [si,~] = rtrcpuls(rollOff, Tau, fs, span);
    symbolsBarker = constBPSK(symbBarker);
    pulseBarker = conv(upsample(symbolsBarker, sps), si);
    
    % Declarations
    foundBarker = 0;
    threshold = 20;
    
    % Start Recording
    record(recording);
    tic;
    win = 1;
    while (true)
        % Record 1 second (entire packet is 0.99s)
        pause(2);
            
        % Fetch Data
        inputSound = getaudiodata(recording, 'single');
        
        t = ((1:length(inputSound))/fs).';
        inputData = inputSound.*(exp(-1i*2*pi*fc*t));
        
        inputData = conv(inputData, lpf);
        
        yt = conv(si, inputData);
        
        if (foundBarker == 1)
            break
        end

        corr = conv(yt, fliplr(pulseBarker));
        
        % Check for correlation with barker code
        [peakV, idxPeak] = max(abs(corr));
        corr(idxPeak-(floor(nBarker/2))*sps:idxPeak+(floor(nBarker/2))*sps) = [];
        
        % Check for difference between barker code and packet
        [peakV2, idxPeak2] = max(abs(corr));
        diff = peakV - peakV2;

        if (peakV > threshold) && (diff > threshold)
            idxPreamble = idxPeak;
            foundBarker = 1;
        end
            
        % If the preamble has not been found in tout sec, return
        if (tout < toc)
            disp('I Found Nothin')
            Xhat=[]; psd=[]; const=[]; eyed=[];            
            return
        end
    end
    stop(recording)
    disp('MF Tic to the TOC')

    % Calculate the delay so that the start of the packet is known    
    delayHat = idxPreamble - length(pulseBarker) + sps*span + 1;
    
    % Remove all bits before the matching in convolution 
    uncutPacket = yt(delayHat:end);
    uncutPacket = uncutPacket(1:(Ns*m+nBarker)*sps);
    
    %% Decision: correct for QPSK and 8PSK
    downPacket = downsample(uncutPacket, sps);     % Downsampling
%     constPilot = downPacket(nBarker+1:nBarker+nPilots);          % Take the first pilot bits, starting at 2 to account for downsampling error
    constBarker = downPacket(1:nBarker);
    constPack = downPacket(nBarker+1:nBarker+Ns);    % The rest of the packet
    
%     stem(real(downPacket));
    
    % Calculating Phaseshift 
%     phaseShift = mean(wrapTo2Pi(angle(constPilot))) - (pi/4);
    
    % Calculating Phaseshift: Barker way
       
    tmpOnes = constBarker(barker == 1);
    tmpOnes = [mean(real(tmpOnes)) mean(imag(tmpOnes))];
    tmpOnes = tmpOnes(:,1) + tmpOnes(:,2)*1i;
    phaseShift = wrapTo2Pi(angle(tmpOnes))-5*pi/4;
    
    %% Matching samples
    % Getting angles from received signal constellation
    constAngle = angle(constPack) - phaseShift;
    % @TODO: Somehow remove for-loop for more efficiency?
    indexSymb = zeros(length(constPack),1); 
    for i = 1:length(constPack)
        % Minimun distance on the constellation
        [~,x] = min(abs(constAngle(i)-QPSKAngle));
        % Matching each symbol with correct constellation
        indexSymb(i) = x;               
    end;
    symbolsRec = indexSymb-1;
    
    %% Bits matched
    bitsGroup = de2bi(symbolsRec, m, 'left-msb');
        
    %% Outputs
    Xhat = reshape(bitsGroup.',[1,m*length(bitsGroup)]);
    const = constPack*exp(-1i*phaseShift);

    tempConst = zeros(1,length(const));
    maxNorm = 0; 
        
    for j=1:length(const)
        if(norm(const(j))>maxNorm)
            maxNorm = norm(const(j));
        end
    end

    for i=1:length(tempConst)    
     tempConst(i)=const(i)./maxNorm;
    end
	const = tempConst;
    
    %% PSD output (WE HAVE TO USE FFT)!!!!!!!!!!! ARRANGE THIS MADAFAKAAAAAAAAAAAAAA
%     XhatdB = 20*log10(abs(const));
%     [psdXhat, fXhat] = pwelch(XhatdB,hamming(128),[],[],fs,'twosided');
%     psd = struct('p',psdXhat,'f',fXhat);
    
    PSDN = length(const);
    fXhat = zeros(1:length(PSDN));
    xdft = fft(const);
    psdx = (1/(2*pi*PSDN))*abs(xdft).^2;
    psdXhat = 20*log10(abs(psdx));
    for i = 1:(PSDN/2-1)
    fXhat(i) = (i*fs)/PSDN;
    end
    for i = PSDN/2 : PSDN
        fXhat(i) =(i*fs)/PSDN - fs;
    end
    psd = struct('p',psdXhat,'f',fXhat);
    %% Eyed output
    %eyed = struct('r',uncutPacket((2+nPilots)*sps:(nPilots+Ns+1)*sps),'fsfd',sps); 
    eyed = struct('r',uncutPacket(sps+nBarker*sps:nBarker*sps+Ns*sps-sps)*exp(-1i*phaseShift),'fsfd',sps);

    %% DEBUGGING ZONE
    % When correcting packet, we loose energy (see Figure 123, the black plot!): ask Keerthi why, he knows
    correctedPacket = uncutPacket(sps+nBarker*sps:nBarker*sps+Ns*sps-sps)*exp(-1i*phaseShift);
    eyediagram(correctedPacket,sps);
    
    scatterplot(const) 
    
    %%
    figure(123)
    subplot(2,1,1)
 
    plot(real(uncutPacket(1+nBarker*sps:nBarker*sps+Ns*sps-sps)))
    title('real')
    hold on
    plot(upsample(real(constPack),sps),'r*')
    plot(upsample(real(const),sps),'g*')
    plot(real(correctedPacket),'k');
    subplot(2,1,2)
  
    plot(imag(uncutPacket(1+nBarker*sps:nBarker*sps+Ns*sps-sps)))
    title('imag')
    hold on
    plot(upsample(imag(constPack),sps),'r*')
    plot(upsample(imag(const),sps),'g*')
    plot(imag(correctedPacket),'k');
    
    %%
%     figure(1);
%     title('Correlation');                       % @TODO: Solve convolution problems so that
%     plot(abs(corr));                            %        we can tun even if sending random bits
%     % Plotting Constellations from received signal
%     scatterplot(constPack*exp(-1i*phaseShift));
%     scatterplot(constPilot*exp(-1i*phaseShift));
% 
%     figure(4); subplot(2,1,1); plot(real(barkerPass), 'b');                         
%                          title('real')
%          subplot(2,1,2); plot(imag(barkerPass), 'r');                        
%                          title('imag')
%                                                
%     figure(23); subplot(2,1,1); plot(real(pulseBarker), 'b');                         
%                          title('real')
%          subplot(2,1,2); plot(imag(pulseBarker), 'r');                        
%                          title('imag')

%     figure(111);
%     subplot(2,1,1);
%     plot(real(yt));
%     subplot(2,1,2);
%     plot(imag(yt));
    %
%end