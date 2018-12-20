%% Network Configuration File

%% Part 1: Basic Setup
% simT: total simulation time 
% Tpkt: time needed for 1 transmission
% Tcont: contention time
% Tupdate: period of piggyback relation update for Original-CHAIN
simT = 10000;
Tpkt = 0.16;
Tdummy = 0.08;
Tcont = 0.1;
Tupdate = 1000;
number_of_mini_slots_per_unit_time = 100; % for arrival process


%% Part 3: Node configuration
N_CHAIN_group = 1;
CHAINs = cell(N_CHAIN_group,1);
CHAINs{1} = [1 2];
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
% For QLB-CSMA
qlen_update_interval = 5;

%% Part 5: Arrival rates for Qth-based algorithm
%rho = 0.97;
arrival_rate = rho*ones(N_nonsat, 1);
%arrival_rate = rho*[15/100, 15/100, 15/100, 10/100, 5/100, 5/100, 5/100, 5/100, 5/100, 5/100];
%arrival_rate_group1 = rho*[100; 90; 80; 70; 60; 50; 40; 30; 20; 10]/100;
%arrival_rate_group2 = rho*10/100*ones(N_half, 1);
%arrival_rate = [arrival_rate_group1; arrival_rate_group2];
%arrival_rate = rho*[50/100, 50/100, 40/100, 40/100, 40/100, 40/100, 40/100, 30/100, 30/100, 30/100, ...
                   % 30/100, 10/100, 10/100, 10/100, 10/100, 10/100, 10/100, 10/100, 5/100, 5/100];
