%% Network Configuration File

%% Part 1: Basic Setup
% simT: total simulation time 
% Tpkt: time needed for 1 transmission
% Tcont: contention time
% Tupdate: period of piggyback relation update for Original-CHAIN
simT = 100000;
Tpkt = 0.2;
Tdummy = 0;
Tcont = 0.5;
Tupdate = 2000;

%% Part 2: Choose algorithm
% Qth-based: 
% Qth-cross-piggyback:
% Qth-plus-Contention:
% Original-CHAIN:
% Original-CHAIN-two-chains:
Mode = 'Original-CHAIN';

%% Part 3: Node configuration
N_CHAIN_group = 1;
CHAINs = cell(N_CHAIN_group,1);
CHAINs{1} = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
DCF = [];
len_CHAIN = cellfun('length',CHAINs);
N_CHAIN_node = sum(len_CHAIN);
N_DCF = length(DCF);
N = N_CHAIN_node + N_DCF;
is_CHAIN = [ones(N_CHAIN_node,1); zeros(N_DCF,1)];
is_DCF = [zeros(N_CHAIN_node,1); ones(N_DCF,1)];
CW_Min = 16;
CW_Max = 1024;
slot_time = 0.009;
DIFS = 0.034;

%% Part 4: Policy parameters
q_threshold = 10; % threshold for reactivation
Tput_shortterm_delta = 0.1;

%% Part 5: Arrival rates for Qth-based algorithm
rho = 0.83;

%arrival_rate = rho*[15/100, 15/100, 15/100, 10/100, 5/100, 5/100, 5/100, 5/100, 5/100, 5/100];
%arrival_rate = rho*25/100*ones(N,1);
arrival_rate = rho*[50/100, 50/100, 40/100, 40/100, 40/100, 40/100, 40/100, 30/100, 30/100, 30/100, ...
                    30/100, 10/100, 10/100, 10/100, 10/100, 10/100, 10/100, 10/100, 5/100, 5/100];
