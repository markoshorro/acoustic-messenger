%% PARAMETERS
% Just the collection of parameters in common for Rx and Tx
%%

%% Constellation
const = [1+1i; -1+1i; -1-1i; 1-1i];

%% Sampling freq
f_samp = 44*10^3;
T_samp = 1/f_samp;

%% Bit rate
R_b = 444;

%% Number of bits to transmit
N = 432;

%% QSPK, 4 different options
m = 2;

%% Alphabet
M = 4; 