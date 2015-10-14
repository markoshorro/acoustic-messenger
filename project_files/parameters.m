%% PARAMETERS
% Just the collection of parameters in common for Rx and Tx

%% Preamble

% GUARD
guard = zeros(1,4);
symbGuard = guard + 1;
nGuard = length(guard);


% Barker code
barker = [1 1 1 1 1 0 0 1 1 0 1 0 1];
symbBarker = barker + 1;
nBarker = length(barker);
% Pilot
pilots = [];
symbPilots = pilots + 1;
nPilots = 10;

%% Constellation
constQPSK = [1+1i; 1-1i; -1+1i; -1-1i]/sqrt(2);
constBPSK = [1+1i -1-1i]/sqrt(2);

QPSKAngle = angle(constQPSK).';
BPSKAngle = angle(constBPSK).';

%% Sampling freq
fs = 44e3;
Ts = 1/fs;

%% Bit rate
Rb = 440;

%% Number of bits to transmit
N = 432;

%% QSPK, 4 different options
m = 2;

%% Alphabet
M = m^2; 

%% Symbol rate
Rs = Rb/m;

%% Symbol time
Tau=1/Rs;

%% Span
span=12;

%% Samples per symbol
sps=round(fs/Rs);

%% Number of symbols per packet
Ns = N/m;

%% Rolloff factor
rollOff = .3;

%% Bwidth: since it is a raised cosine
bandWd = (1+rollOff)/(2*Tau);
Wn = 2*bandWd/fs;
cutOff = 150;
