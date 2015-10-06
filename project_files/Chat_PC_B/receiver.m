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
    run('../parameters.m')

    %% Audio data collection
    message = zeros(1,1000) + 0.5;          %testing dummy
    recording = audiorecorder(44000,8,1);   %Creating recording Object
    record(recording);                      %start recording
    tic;        %start counter, to keep track of recording time of each
                %recording segment
    
    while message(end) == 0.5;  %marker condition (dummy)
        while toc < 2;          %waits for 'toc' seconds to record     
        end                 
        pause(recording);       %Pause recording
        message = getaudiodata(recording,'single'); %fetch data
        resume(recording);                         
    end    
    stop(recording);    %stop recording after finding correct packet size
    
    %% Passband to baseband
    t = (0:1/length(message):1-1/length(message)).';

    data = message; %%%%%%%%%%%%%%%% just for testing
    data = data.*(exp(-1i*2*pi*fc*t));
    
    %% Demodulation (MF)
    [si,~] = rtrcpuls(0.3, Tau, fs, span);
    yt = conv(si, data);
    yt = yt(sps*span:end-sps*span);
    
    %% Decision: correct for QPSK and 8PSK
    const = downsample(yt, sps);
    yangle = angle(const);
    constangle = angle(constQPSK).';
    index_symb=zeros(length(const),1);
    for i=1:length(const)
        [~,x] = min(abs(yangle(i)-constangle));
        index_symb(i)=x;
    end;
    symbols_rec = index_symb-1;
    
    bits_group = de2bi(symbols_rec);
    %% XHAT output
    Xhat = reshape(bits_group.',[1,m*length(bits_group)]);
    
    Xhat
    const
    figure(14);
    subplot(2,1,1);
    plot(real(yt));
    subplot(2,1,2);
    plot(imag(yt));
    
    
    %% PSD output
    Xhat_dB = 20*log10(Xhat);
    [psd_Xhat, f_Xhat] = pwelch(Xhat_dB,fs); %psd_Xhat needs to be normalized so the max reaches 0 dB
    psd_Xhat
    f_Xhat
    field1 = 'p';
    field2 = 'f';
    psd = struct(field1,psd_Xhat,field2,f_Xhat);
    
    %% eyed output
    field3 = 'fsfd';
    field4 = 'r';
    eyed = struct(field4,yt,field3,sps);
    
end