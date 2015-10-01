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
    fc = 5000;
   
    
    %% Audio data collection 
    
    message = zeros(1,1000) + 0.5;          %testing dummy
    recording = audiorecorder(44000,8,1);   %Creating recording Object
    record(recording);                      %start recording
    tic;        %start counter, to keep track of recording time of each
                %recording segment
    
    while message(end) == 0.5;  %marker condition (dummy)
        while toc < 1;          %waits for 'toc' seconds to record     
        end                 
        pause(recording);       %Pause recording
        message = getaudiodata(recording,'single'); %fetch data
        resume(recording);                         
    end    
    stop(recording);    %stop recording after finding correct packet size
    
    %message = dec2bin(message)';
    %message = reshape(message,1,8*length(message));
    
    Xhat = message;
    %scatterplot(Xhat)    
    
    %% Passband to baseband
    t = (0:1/length(Xhat):1-1/length(Xhat)).';
    
    %data = s_passband; % just for testing
    data = Xhat.*(exp(-1i*2*pi*fc*t)/sqrt(2));
    
    %% Demodulation (MF)
    [si,~] = rtrcpuls(0.3, Tau, fs, span);
    Y = conv(si, data);
    
    figure(1); plot(real(Y))
    figure(2); plot(downsample(real(Y),200))
    scatterplot(downsample(Y,200))
    
    %% Symbol to bits
    % we should have a 1D vector with values between [1,4]
    
    %symbols = buffer(message, m)';
    %symbols = bin2dec(symbols);
    %Xhat = constQPSK(symbols + 1);
    %scatterplot(Xhat)
%end