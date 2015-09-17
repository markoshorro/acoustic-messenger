function transmitter(packet,fc)

f_samp = 44*10^3; %sampling freq
T_samp= 1/f_samp;
R_b = 444; %bit rate


N = 432; % number of bits to transmit

m= 2; %QSPK, 4 different options
M=4; %Alphabet

R_s = R_b/m; %symbol rate
f_srs = f_samp/R_s; %number of samples in each symbol

Tau=1/R_s;



QPSK_const = [1+1i; -1+1i; -1-1i; 1-1i];


bits_group = buffer(packet,m)'; % dela upp 2 och 2 (reshape funkar också)

messages = bi2de(bits_group)+1; %kalla varje kombination en siffra

symbols = QPSK_const(messages) %anknyt varje komb med vår konstellation


scatterplot(symbols);

symbols_up = upsample(symbols,f_srs); 

[y,t]=rtrcpuls(0.3,Tau,f_samp,3);



A=conv(y,symbols_up); 
                                        


end