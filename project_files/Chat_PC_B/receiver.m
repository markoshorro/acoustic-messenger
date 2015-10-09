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
    
    % Testing
    close all;
    fc = 6000;
    
    %% Audio data collection
    channels = 1;       % Mono
    recordBits = 16;    % More precision
    
    % Creating recording object
    recording = audiorecorder(fs, recordBits, channels);   
    
    %% Low pass filter params
    order = 2^10;
    lpf = fir1(order, Wn);   

    %% @TODO
    % RECORDING AUDIO!
    % Create the barker pulse train
    [si,~] = rtrcpuls(rollOff, Tau, fs, span);
    symbolsBarker = constBPSK(symbBarker);
    pulseBarker = conv(upsample(symbolsBarker, sps), si);
    
    % Start Recording
    record(recording);
    
%     found = false;
%     recording = audiorecorder(fs, recordBits, channels);   % Creating recording Object
%     record(recording);              % start recording
%     index = 1;
%     win = nBarker*sps;
%     tic;
%     while (found==false)
%         % Recording a window
%         pause(win/fs);
%
%         % Getting data from mic
%         in = getaudiodata(recording, 'single');
%
%         % Select a window
%         tmp = in(index:index+win);
%
%         % Update index:
%         index = index + win;
%
%         % Finding the preamble
%         baseIndex = conv(si,index.*(exp(-1i*2*pi*fc*t));
%         corr = conv(baseIndex, fliplr(pulseBarker));
%         [peak, indexPeak] = max(corr);         
%
%         if (peak > threshold)
%               indexPreamble = indexPeak; 
%               found = true;
%         end
%         
%         % If the preamble has not been found in tout sec, return
%         if (tout<toc)
%             Xhat=[]; psd=[]; const=[]; eyed=[];
%             return;
%         end     
%     end
%    pause(1);
%    inputSound = getaudiodata(recording,'single');
%    stop(recording);
%     
    pause(4);
    inputSound = getaudiodata(recording,'single');
    
    % Stop recording after finding correct packet size
    stop(recording);                
    
    %% Passband to baseband
    t = ((1:length(inputSound))/fs).';
    inputData = inputSound.*(exp(-1i*2*pi*fc*t));
    
    %% Low pass filter 
    inputData = conv(inputData, lpf);
    
    %% Demodulation (MF)
    % Convolution between received signal and rtrc pulse
    yt = conv(si, inputData);
    
    %% TO ARRANGE. NOT WORKING PROPERLY
    % Autocorrelate signal with barker code
    corr = conv(yt, fliplr(pulseBarker));
    
    % Finding autocorrelation peak  
    [~,startSamples] = max(corr);  
    
    %startSamples = indexPreamble + nPilots*sps + nBarker*sps;

    % Calculate the delay so that the start of the packet is known    
    delay_hat = startSamples - nBarker*sps;
    
    % Remove all bits before the matching in convolution 
    uncutPacket = yt(delay_hat:end);
    
    %% Decision: correct for QPSK and 8PSK
    downPacket = downsample(uncutPacket, sps);     % Downsampling
    constPilot = downPacket(2:nPilots+1);          % Take the first pilot bits, starting at 2 to account for downsampling error
    constPack = downPacket(nPilots+2:nPilots+Ns+1);    % The rest of the packet
    
    % Calculating Phaseshift
    phaseShift = wrapTo2Pi(mean(angle(constPilot))) - (pi/4);
    
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
    bitsGroup = de2bi(symbolsRec, m);
        
    %% XHAT output
    Xhat = reshape(bitsGroup.',[1,m*length(bitsGroup)]);
    
    %% PSD output
    XhatdB = 20*log10(constPack);
    [psdXhat, fXhat] = pwelch(abs(XhatdB),...
                                hamming(128),[],[],fs,'twosided');
    psd = struct('p',psdXhat,'f',fXhat);
    
    %% Eyed output
    eyed = struct('r',yt,'fsfd',sps);
    
    %% DEBUGGING ZONE
    figure(1);
    title('Correlation');                       % @TODO: Solve convolution problems so that
    plot(abs(corr));                            %        we can tun even if sending random bits
    % Plotting Constellations from received signal
    scatterplot(constPack*exp(-1i*phaseShift));
    scatterplot(constPilot*exp(-1i*phaseShift));

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
    
%end