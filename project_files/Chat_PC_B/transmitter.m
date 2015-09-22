%load '../parameters.m'
span=6;
R_b=440;
R_s=R_b/m;
f_srs=fs/R_s;
packet = randsrc(1,N,[0 1]);        % Just for test

bits_group = buffer(packet,m)';     % split 2 and 2 (the function reshape also works)
messages = bi2de(bits_group)+1;     % call each combination a number
symbols = const(messages);          % match each number with our constellation

symbols_up = upsample(symbols,round(f_srs)); % Space the symbols fsfd apart, to enable pulse shaping using conv.
[si,t] = rtrcpuls(0.3, Tau, f_samp, span);
s=conv(si,symbols_up); 
                    
% figure; subplot(2,1,1); plot(real(s), 'b');                         
%                         title('real')
%         subplot(2,1,2); plot(imag(s), 'b');                        
%                         title('imag')
s_tailness = s(f_srs*span:end-f_srs*span);

% figure; subplot(2,1,1); plot(real(s_tailness), 'b');                         
%                         title('real')
%         subplot(2,1,2); plot(imag(s_tailness), 'b');                        
%                         title('imag')

 t=(0:1/fs:((length(symbols)/216)-1/fs)).';
 t = t(1:length(s_tailness));

 fc= 5000;                       
t_tx = 1:length(s_tailness);
s_passband = real(s_tailness.*(sqrt(2)*exp(1i*2*pi*fc*t)));
figure;
pwelch(s_passband,hamming(512),[],[],f_samp,'centered');