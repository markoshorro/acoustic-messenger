function transmitter(packet, fc)
	%% TRANSMITTER FUNCTION
    % Group 13
    % Introduction to Communication Engineering. September 2015 
    %
    % Implementation of the receiver.
    % INPUT:  packet  - vector with data
    %         fc    - carrier frequency
    %
    % OUTPUT: SOUND to the speaker (
    %
    
    run('../parameters.m');
    
%     packet = randsrc(1,N,[0 1]);        % Just for test

    bits_group = buffer(packet,m)';     % split 2 and 2 (the function reshape also works)
    messages = bi2de(bits_group)+1;     % call each combination a number
    symbols = constQPSK(messages);          % match each number with our constellation

    symbols_up = upsample(symbols,round(sps)); % Space the symbols fsfd apart, to enable pulse shaping using conv.
    [si,~] = rtrcpuls(0.3, Tau, fs, span);
    s=conv(si,symbols_up); 
    
%     figure; subplot(2,1,1); plot(real(s), 'b');                         
%                              title('real')
%              subplot(2,1,2); plot(imag(s), 'b');                        
%                              title('imag')
    s_tailness = s(sps*span:end-sps*span);

%     figure; subplot(2,1,1); plot(real(s_tailness), 'b');                         
%                              title('real')
%              subplot(2,1,2); plot(imag(s_tailness), 'b');                        
%                              title('imag')

    t = (0:1/fs:((length(symbols)/216)-1/fs)).';
    t = t(1:length(s_tailness));
                      
    s_passband = real(s_tailness.*(sqrt(2)*exp(1i*2*pi*fc*t)));
%     figure(1);
%     pwelch(s_passband,hamming(512),[],[],fs,'centered');
%     
%     figure(2)
%     plot(s_passband)
    s_passband = s_passband/max(s_passband);
    sound(s_passband,44e3);
  
end