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
    %close all;
    % Testing
    fc = 6000;
    
    %% Audio data collection
    channels = 1;   
    recordBits = 16;
    
    % Preamble stuff
%     [si,~] = rtrcpuls(rollOff, Tau, fs, span);
%     symbolsBarker = constBPSK(symbBarker);
%     pulseBarker = conv(upsample(symbolsBarker, round(sps)), si);

%     figure(4); subplot(2,1,1); plot(real(barkerPass), 'b');                         
%                          title('real')
%          subplot(2,1,2); plot(imag(barkerPass), 'r');                        
%                          title('imag')
%                                                
%     figure(23); subplot(2,1,1); plot(real(pulseBarker), 'b');                         
%                          title('real')
%          subplot(2,1,2); plot(imag(pulseBarker), 'r');                        
%                          title('imag')
    
%     message = zeros(1,1000) + 0.5;  % testing dummy
    recording = audiorecorder(fs, recordBits, channels);   % Creating recording Object
    record(recording);              % start recording
    tic;                            % start counter, to keep track of recording 
                                    %  time of each recording segment
    %% @TODO
    % RECORDING AUDIO! 
%     found = false;
%     while (found==false)
%         % Getting data from mic
%         r = getaudiodata(recording, 'single');
%
%         % Finding the preamble @TODO
%         
%         % If the preamble has not been found in tout sec, return
%         if (tout<toc)
%             Xhat=[]; psd=[]; const=[]; eyed=[];
%             return;
%         end     
%     end
        
      %pause(1);      
%     while message(end) == 0.5;    % marker condition (dummy)
        while toc < 4;              % waits for 'toc' seconds to record     
        end
        
        
        
%        pause(recording);          % Pause recording
%        message = getaudiodata(recording,'single'); %fetch data
%        resume(recording);                         
%     end
    message = getaudiodata(recording,'single');
    
%     N = max(1024,length(message));
%     P = fftshift(fft(message,N));     
%     fvec = (fs/N).*[0:(N-1)];
%     fvec = (fvec < fs/2); 
%    % figure(8); plot(fvec,20*log10(abs(P)));  
%     
%     xdftYouCareAbout = abs(P(1:round(N/2))); % Take the absolute magnitude.
% 
% [maxVal, index] = max(xdftYouCareAbout); % maxVal is your (un-normalized) maximum amplitude
% 
% maxFreq = freqsYouCareAbout(index); % This is the frequency of your dominant signal. 
    
%     tesdf= 20*log10(abs(P))
%    
%     [~,b] = max(fvec);
%     b

    stop(recording);                % Stop recording after finding correct packet size
    
    %% Passband to baseband
    t = ((1:length(message))/fs).';
    data = message.*(exp(-1i*2*pi*fc*t));
    
%     pwelch(data, hamming(128),[],[],fs,'twosided');
    
    %% Low pass filter
    order = 2^10;
    Wn = 2*cutOff/fs;
    lpf = fir1(order, Wn);    
    %data = conv(data, lpf);
    
%     pwelch(data, hamming(128),[],[],fs,'twosided');
    
    %% Demodulation (MF)
    yt = conv(si, data);            % Convolution between received signal and rtrc pulse
    %yt = yt(sps*span:end-sps*span);

    [si,~] = rtrcpuls(rollOff, Tau, fs, span);
    symbolsBarker = constBPSK(symbBarker);
    const = downsample(yt, sps);        % Downsampling
    pulseBarker = conv(symbolsBarker, si);
    
    % TO TRY -> t   his should get the peak value
    corr = conv(const, fliplr(symbolsBarker));
    figure(122);
    title('Correlation');
    plot(abs(corr));
    
        [~,PackStart] = max(abs(corr));

        PackStart = PackStart
        
    %     figure(111);
    %     subplot(2,1,1);
    %     plot(real(yt));
    %     subplot(2,1,2);
    %     plot(imag(yt));

        %% Decision: correct for QPSK and 8PSK

        constPack = const(PackStart+10:PackStart+10+Ns-1);
        constPilot = const(PackStart:PackStart+9);
        constBarker = const(PackStart-nBarker:PackStart-1);

        scatterplot(constPack);                 % Plotting Constellations from received signal
        scatterplot(constBarker);
        scatterplot(constPilot);
        constAngle = angle(constPack);          % Getting angles from received signal constellation
        indexSymb = zeros(length(constPack),1); % Declaration before for-loop

        for i = 1:length(constPack)             % @TODO: Somehow remove for-loop for more efficiency?
            [~,x] = min(abs(constAngle(i)-QPSKAngle)); 
            indexSymb(i) = x;               % Matching each symbol with correct constellation
        end;

        symbolsRec = indexSymb-1;

        bitsGroup = de2bi(symbolsRec);


        h = 1:length(yt)/length(constPack):length(yt);
%     figure(123)
%     plot(real(yt),'k')
%     hold on
%     plot(h,real(constPack),'*r')
%     hold off
    
    %% XHAT output
    Xhat = reshape(bitsGroup.',[1,m*length(bitsGroup)]);
    
    %% PSD output
    XhatdB = 20*log10(constPack);
    [psdXhat, fXhat] = pwelch(abs(XhatdB),...
                                hamming(128),[],[],fs,'twosided');
    psd = struct('p',psdXhat,'f',fXhat);
    
    %% eyed output
    eyed = struct('r',yt,'fsfd',sps);
%end