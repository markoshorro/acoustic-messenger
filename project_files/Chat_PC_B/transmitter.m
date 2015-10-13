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
    global symbols;
    run('../parameters.m');
   % pilot = zeros(1,20);
    packet = [packet'];        % Just for test
%     fc = 6000;

    % Split in m columns
    bitsGroup = buffer(packet,m)';     
    messages = bi2de(bitsGroup,'left-msb')+1;
    
    % Match each number with our constellation
    symbols = constQPSK(messages);
    
    % Match barker code with BPSK constellation
    symbolsBarker = constBPSK(symbBarker);
    symbols = [symbolsBarker.'; symbols];
    
    % Space the symbols fsfd apart, to enable pulse shaping using conv
    symbolsUp = upsample(symbols, round(sps));
        
    % Pulse convolution
    [si,~] = rtrcpuls(rollOff, Tau, fs, span);
    st = conv(si, symbolsUp); 
    
    % Getting convoluted signal without heads nor tails
    sTailless = st(sps*span:end-sps*span);
    
    % Converting onto passband signal
    t = ((1:length(sTailless))/fs).';
    sPassband = real(sTailless.*(exp(1i*2*pi*fc*t)));

    % Normalized values
    sPassband = sPassband/max(sPassband);
    
    % Output through the speaker
    sound(sPassband, fs);

    %% DEBUGGING
%     figure(1); subplot(2,1,1); plot(real(sTailless), 'b');                         
%                              title('Real Part Concatinated - BaseBand')
%              subplot(2,1,2); plot(imag(sTailless), 'r');                        
%                              title('Imaginary Part Concatinated - BaseBand')
                             
%     figure(2); subplot(2,1,1); plot(real(st), 'b');                         
%                              title('Real Part - BaseBand')
%              subplot(2,1,2); plot(imag(st), 'r');                        
%                              title('Imaginary Part - BaseBand')

%     figure(3); subplot(2,1,1); plot(real(sPassband), 'b');                         
%                              title('Real Part - PassBand')
%              subplot(2,1,2); plot(imag(sPassband), 'r');                        
%                              title('Imaginary Part PassBand')
%     figure(1);
%     pwelch(sPassband,hamming(512),[],[],fs,'centered');
%     
%     figure(2)
%     plot(sPassband)
%     figure; subplot(2,1,1); plot(real(sTailness), 'b');                         
%                              title('real')
%              subplot(2,1,2); plot(imag(sTailness), 'b');                        
%                              title('imag')
end