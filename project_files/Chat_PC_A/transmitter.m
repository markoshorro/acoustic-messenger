%function transmitter(packet,fc)
packet = randsrc(1,N,[0 1]);
fc=8000

f_samp = 44*10^3; %sampling freq
T_samp= 1/f_samp;
R_b = 444; %bit rate

disp('New test')

N = 432; % number of bits to transmit

m= 2; %QPSK, 4 different options
M=4; %Alphabet

R_s = R_b/m; %symbol rate
f_srs = f_samp/R_s; %number of samples in each symbol

Tau=1/R_s;



QPSK_const = [1+1i; -1+1i; 1-1i; -1-1i]/sqrt(2);


bits_group = buffer(packet,m)'; % split 2 and 2 (the function reshape also works)

messages = bi2de(bits_group)+1; %call each combination a number

symbols = QPSK_const(messages); %match each number with our constellation


%scatterplot(symbols);

symbols_up = upsample(symbols,round(f_srs)); 

[y,t]=rtrcpuls(0.3,Tau,f_samp,3);



A=conv(y,symbols_up); 
                    
figure; subplot(2,1,1); plot(real(A), 'b');                         
                        title('real')
        subplot(2,1,2); plot(imag(A), 'b');                        
                        title('imag')
                        
                        



A= [real(A).*cos(2*pi*fc) imag(A).*sin(2*pi*fc)];
size(A)



N=length((real(A))); %length of one vector, i.e our window length (time)
%N=pow2(nextpow2(L)); %length of our fft (total bins)
N=1024;

f_vec=(f_samp/N)*(-floor(N/2):1:ceil(N/2)-1);

% figure(89)
% plot(f_vec,fft(A,N))

length(f_vec)
length(fft(A,N))
figure(50)
plot(f_vec,20*log10((abs(fftshift(fft(A,N)))).^2))

sound(real(A),f_samp)





%end