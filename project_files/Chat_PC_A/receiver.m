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
    fc = 4000;
    
    %% Audio data collection
    channels = 1;   
    recordBits = 8;
    
    % Preamble stuff
    [si,~] = rtrcpuls(rollOff, Tau, fs, span);
    symbolsBarker = constBPSK(symbBarker);
    pulseBarker = conv(upsample(symbolsBarker, round(sps)), si);
    t = ((1:length(pulseBarker))/fs);
    barkerPass = real(pulseBarker.*(exp(1i*2*pi*fc*t)));
    barkerPass = barkerPass/max(barkerPass);

    figure(4); subplot(2,1,1); plot(real(barkerPass), 'b');                         
                         title('real')
         subplot(2,1,2); plot(imag(barkerPass), 'r');                        
                         title('imag')
    
    message = zeros(1,1000) + 0.5;          %testing dummy
    recording = audiorecorder(fs, recordBits, channels);   %Creating recording Object
    record(recording);                      %start recording
    tic;        %start counter, to keep track of recording time of each
                %recording segment
%     found = false;
%     while (found==false)
%         % We have to find the preamble
%         r = getaudiodata(recording, 'single');
%         
%         % 
%         
%         % If the preamble has not been found, return
%         if (tout<toc)
%             Xhat=[]; psd=[]; const=[]; eyed=[];
%             return;
%         end     
%     end
        
                
%     while message(end) == 0.5;  %marker condition (dummy)
        while toc < 1;          %waits for 'toc' seconds to record     
        end                 
%        pause(recording);       %Pause recording
%        message = getaudiodata(recording,'single'); %fetch data
%        resume(recording);                         
%     end
    message = getaudiodata(recording,'single');
    stop(recording);    %stop recording after finding correct packet size
    
    corr = conv(message, fliplr(barkerPass));
    plot(corr);
    
    %% Passband to baseband
    t = ((1:length(message))/fs).';
    data = message.*(exp(-1i*2*pi*fc*t));
    
    %% Demodulation (MF)
    yt = conv(si, data);
    yt = yt(sps*span:end-sps*span);
    
    %% Decision: correct for QPSK and 8PSK
    const = downsample(yt, sps);
    yangle = angle(const);
    constAngle = angle(constQPSK).';
    indexSymb=zeros(length(const),1);
    for i=1:length(const)
        [~,x] = min(abs(yangle(i)-constAngle));
        indexSymb(i)=x;
    end;
    symbolsRec = indexSymb-1;
    
    bitsGroup = de2bi(symbolsRec);
    
    %% XHAT output
    Xhat = reshape(bitsGroup.',[1,m*length(bitsGroup)]);
    
    %% PSD output
    XhatdB = 20*log10(const);
    [psdXhat, fXhat] = pwelch(abs(XhatdB),...
                                hamming(128),[],[],fs,'twosided');
    psd = struct('p',psdXhat,'f',fXhat);
    
    %% eyed output
    eyed = struct('r',yt,'fsfd',sps);
%end