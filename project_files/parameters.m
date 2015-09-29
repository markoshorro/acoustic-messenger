%% PARAMETERS
% Just the collection of parameters in common for Rx and Tx
%%

%% Constellation
constQPSK = [1+1i; -1+1i; -1-1i; 1-1i];

%% Sampling freq
fs = 44*10^3;
Ts = 1/fs;

%% Bit rate
Rb = 440;

%% Number of bits to transmit
N = 432;

%% QSPK, 4 different options
m = 2;

%% Alphabet
M = 4; 

%% Symbol rate
Rs = Rb/m;

%% Symbol time
Tau=1/Rs;

%% Span
span=6;

%% Samples per symbol
sps=fs/Rs;
