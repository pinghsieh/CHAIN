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
% Qth-cross-piggyback:
% Qth-plus-Contention:

Mode = 'Qth-plus-Contention';

%% Part 3: Node configuration
N_CHAIN_group = 1;
CHAINs = cell(N_CHAIN_group,1);
CHAINs{1} = [1 2];
DCF = 3;
len_CHAIN = cellfun('length',CHAINs);
N_CHAIN_node = sum(len_CHAIN);
N_DCF = length(DCF);
N = N_CHAIN_node + N_DCF;
is_CHAIN = [ones(N_CHAIN_node,1); zeros(N_DCF,1)];
is_DCF = [zeros(N_CHAIN_node,1); ones(N_DCF,1)];

%% Part 4: Policy parameters
q_threshold = 100; % threshold for reactivation

%% Part 5: Arrival rates for Qth-based algorithm
rho = 0.95;

arrival_rate = rho*[1/5, 1/5, 1/5];


