%% Network Configuration File

%% Part 1: Basic Setup
% simT: total simulation time 
% Tpkt: time needed for 1 transmission
% Tcont: contention time
simT = 500000;
Tpkt = 1;
Tdummy = 0;
Tcont = 1;

%% Part 2: Choose algorithm
% Qth-based: 
% Qth-based-DCF-priority:
% Qth-based-Chain-wide-contention:
% Qth-cross-piggyback:
% Qth-plus-Contention:
% Qth-plus-Contention-DCF-priority
% Pure-DCF
Mode = 'Qth-plus-Contention-DCF-priority';

%% Part 3: Node configuration
N_CHAIN_group = 1;
CHAINs = cell(N_CHAIN_group,1);
CHAINs{1} = [1 2 3];
DCF = 4;
len_CHAIN = cellfun('length',CHAINs);
N_CHAIN_node = sum(len_CHAIN);
N_DCF = length(DCF);
N = N_CHAIN_node + N_DCF;
is_CHAIN = [ones(N_CHAIN_node,1); zeros(N_DCF,1)];
is_DCF = [zeros(N_CHAIN_node,1); ones(N_DCF,1)];

%% Part 4: Policy parameters
q_threshold = 10; % threshold for reactivation

%% Part 5: Arrival rates for Qth-based algorithm
rho = 0.99;

arrival_rate = rho*[1/100, 1/100, 1/100, 48/100];


