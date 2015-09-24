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
    load '../parameters.m'
    RecordingTime = 6; %How long (in seconds) we will record sound
    
    %% Audio data collection 
    
    recording = audiorecorder(44000,8,1);
    record(recording);
    T = timer('TimerFcn',@(~,~)disp('Test running'),'StartDelay',RecordingTime);
    start(T)    %waiting for 'RecordingTime' seconds
    wait(T)
    stop(recording);
    data = getaudiodata(recording,'int8');
    
    
    %% Passband to baseband
    % Where Vt
    ak = Yq*sqrt(2)*cos(2*pi*fc*t);
    bk = Yq*sqrt(2)*sin(2*pi*fc*t);
    
    %% Syncronization
    
    %% Demodulation
    
    %% Symbol to bits
    % we should have a 1D vector with values between [1,4]
    bits_group = de2bi(message, m);
    Xhat = buffer(bits_group, 1);
    
    
end