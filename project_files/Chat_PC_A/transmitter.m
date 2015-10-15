function transmitter( packet, fc )
	%% TRANSMITTER FUNCTION
    % Group 13
    % Introduction to Communication Engineering. September 2015 
    %
    % Implementation of the receiver.
    % INPUT:  packet  - vector with data
    %         fc    - carrier frequency
    %
    % OUTPUT: SOUND to the speaker
    %
    
    run('../parameters.m');
    packet = [packet'];        % Just for test

    % Split in m columns
    bitsGroup = buffer(packet,m)';     
    messages = bi2de(bitsGroup,'left-msb')+1;
    
    % Match each number with our constellation
    symbols = constQPSK(messages);
    
    % Match barker code with BPSK constellation
    symbolsBarker = constBPSK(symbBarker);
    
    % Match guard symbol with BPSK constellation
    symbolsGuard = constBPSK(singleGuard); 

    % Concatenate all symbols
    symbols = [symbolsBarker.'; symbolsGuard.'; symbols];
    
    % Space the symbols fsfd apart, to enable pulse shaping using conv
    symbolsUp = upsample(symbols, sps);
        
    % Pulse convolution
    [si,~] = rtrcpuls(rollOff, Tau, fs, span);
    st = conv(si, symbolsUp); 
    
    % Getting convoluted signal without heads nor tails
    sTailless = st(sps*span:end-sps*span);
    sTailless = [sTailless(1:sps*nBarker); upsample(guard,sps).';...
        sTailless(sps*nBarker+1:end)];
    
    % Converting onto passband signal
    t = ((1:length(sTailless))/fs).';
    sPassband = real(sTailless.*(exp(1i*2*pi*fc*t)));

    % Normalized values
    sPassband = sPassband/max(sPassband);
    
    % Output through the speaker
    sound(sPassband, fs);
end
