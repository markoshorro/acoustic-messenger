%% PARAMETERS
% Just the collection of parameters in common for Rx and Tx
%%

%% Constellation
const = [1+1i; -1+1i; -1-1i; 1-1i];

%% Sampling freq
fs = 44*10^3;
T_samp = 1/fs;

%% Bit rate
R_b = 440;

%% Number of bits to transmit
N = 432;

%% QSPK, 4 different options
m = 2;

%% Alphabet
M = 4; 

%% Symbol rate
R_s = R_b/m;

%% Symbol time
Tau=1/R_s;

%% Span
span=6;

%% Samples for symbol
f_srs=fs/R_s;
