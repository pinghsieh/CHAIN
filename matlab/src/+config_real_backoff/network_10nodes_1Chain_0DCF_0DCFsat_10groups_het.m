%% Network Configuration File

%% Part 1: Basic Setup
% simT: total simulation time 
% Tpkt: time needed for 1 transmission
% Tcont: contention time
% Tupdate: period of piggyback relation update for Original-CHAIN
simT = 30000;
Tpkt = 0.16;
Tdummy = 0.08;
Tcont = 0.1;
Tupdate = 1000;
number_of_mini_slots_per_unit_time = 100; % for arrival process


%% Part 3: Node configuration
N_CHAIN_group = 1;
CHAINs = cell(N_CHAIN_group,1);
CHAINs{1} = 1:1:10;
DCF = [];
DCFsat = [];
len_CHAIN = cellfun('length',CHAINs);
N_CHAIN_node = sum(len_CHAIN);
N_DCF = length(DCF);
N_DCFsat = length(DCFsat);
N = N_CHAIN_node + N_DCF + N_DCFsat;
N_nonsat = N_CHAIN_node + N_DCF;
N_half = round(N_nonsat/2);
is_CHAIN = [ones(N_CHAIN_node,1); zeros(N_DCF,1); zeros(N_DCFsat,1)];
is_DCF = [zeros(N_CHAIN_node,1); ones(N_DCF,1); zeros(N_DCFsat,1)];
is_DCFsat = [zeros(N_CHAIN_node,1); zeros(N_DCF,1); ones(N_DCFsat,1)];
CW_Min = 16;
CW_Max = 1024;
slot_time = 0.009;
DIFS = 0.034;
p_star = 1;
channel_pn = p_star*ones(N,1);

%% Part 4: Policy parameters
q_threshold = 10; % threshold for reactivation
Tput_shortterm_delta = 0.4;
dummy_packet_limit = 20;

%% Part 5: Arrival rates for Qth-based algorithm
N_arrival_groups = 10;
arrival_rate = rho*[10; 9; 8; 7; 6; 5; 4; 3; 2; 1];

