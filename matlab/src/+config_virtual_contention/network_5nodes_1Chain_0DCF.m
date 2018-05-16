%% Network Configuration File
MAX_TOL = 1e-10;
%% Part 1: Basic Setup
% simT: total simulation time 
% Tpkt: time needed for 1 transmission
% Tcont: contention time
% Tupdate: period of piggyback relation update for Original-CHAIN
% time unit = 1ms
simT = 100000;
Tpkt = 0.5;
Tdummy = 0;
Tcont = 0.5;
Tupdate = 10;
number_of_mini_slots_per_unit_time = 100;

%% Part 2: Choose algorithm
% Qth-based: 
% Qth-cross-piggyback:
% Qth-plus-Contention:
% Original-CHAIN:
% Original-CHAIN-two-chains:
% Pure-DCF:
Mode = 'Pure-DCF';


%% Part 3: Node configuration
N_CHAIN_group = 1;
CHAINs = cell(N_CHAIN_group,1);
CHAINs{1} = [1 2 3 4 5];
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
Tput_shortterm_delta = 0.4;

%% Part 5: Arrival rates for Qth-based algorithm
rho = 0.5;

%arrival_rate = rho*[80/100, 50/100, 40/100, 20/100, 10/100];
arrival_rate = rho*[40/100, 40/100, 40/100, 40/100, 40/100];

%arrival_rate = rho*[25/100, 20/100, 15/100, 10/100, 5/100];
%arrival_rate = rho*[1/6, 1/6, 1/6, 1/6, 1/6];

