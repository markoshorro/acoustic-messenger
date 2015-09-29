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
    
    recording = audiorecorder(44000,8,1);
    record(recording);
    message = zeros(1,1000) + 0.5;  %testing dummy
    tic;
    
    while message(end) == 0.5; %marker condition
        while toc < 1;         
        end
        pause(recording);
        message = getaudiodata(recording,'uint8');        
        resume(recording);
        
    end    
    stop(recording);  
    
    
    %% Passband to baseband
  
    data = real(data.*(exp(-1i*2*pi*fc*t)/sqrt(2)));
    
    %% LP-filter
    
    [b,a] = butter(2,2*fc/fs,'low'); %fc/(fs/2) -> normalized fc
    
    %h=fvtool(b,a) %if you want to look at it
    
    Y = filter(b,a,data);
    
    
    %% Syncronization
    
    %% Demodulation
    [si,~] = rtrcpuls(0.3, Tau, fs, span);
    %Y=conv(si,Y);
    
    
    Y = filter(si,1,Y)
    
    plot(Y)
    
     %% Symbol to bits
    % we should have a 1D vector with values between [1,4]
    message = dec2bin(message)';
    message = reshape(message,1,8*length(message));
    Symbols = buffer(message, m)';
    
end