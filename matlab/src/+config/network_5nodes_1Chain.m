%% Network Configuration File

%% Part 1: Basic Setup
% simT: total simulation time 
% Tpkt: time needed for 1 transmission
% Tcont: contention time
simT = 100000;
Tpkt = 1;
Tyield = 1;
Tcont = 1;

%% Part 2: Choose algorithm
% Qth-based: 
% Qth-cross-piggyback:
% Qth-plus-Contention:

Mode = 'Qth-based';

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
round_count = zeros(2^N, 1);

%% Part 4: Arrival rates for Qth-based algorithm
arrival_rate = 1*[0.25, 0.2, 0.1, 0.1]/1.035;
%arrival_rate = 0.96*[1/15*0.3, 1/10*0.2, 1/15*0.4, 1/15*0.05, 1/15*0.05];
%arrival_rate = 0.96*[1/15*0.15, 1/10*0.1, 1/15*0.4, 1/15*0.2, 1/15*0.15];
%arrival_rate = 1*[0.04, 0.01, 0.33, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01];
%arrival_rate = 1*[0.04, 0.02, 0.2, 0.2, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01];
%arrival_rate = 0.995*[1/15*0.8, 1/15*0.15, 1/15*0.05];

%arrival_rate = 1*[1/15*0.399, 1/10*0.15, 1/10*0.25, 1/15*0.199, 0];
%arrival_rate = 0.95*1/100*ones(N,1);
%arrival_rate = 1.01*1/55*ones(N,1);
%arrival_rate = 0.99*[1/30, 1/40, 1/50, 1/250, 1/1000];
%arrival_rate = 0.99*[1/50, 1/1000, 1/40, 1/250, 1/30];
%arrival_rate = [1/16, 1/500, 1/600, 1/1000, 1/1000];
%arrival_rate = 0.99*[1/55, 1/55, 1/55, 1/55, 1/55];
%arrival_rate = 1/85*ones(N,1);
%arrival_rate = 0.995*[1/20, 1/50, 1/1000, 1/1000, 1/1000, 1/1000, 1/2000, 1/2000];
%arrival_rate = 0.995*[1/20, 1/100, 1/200, 1/200, 1/500, 1/1000, 1/2000, 1/2000, 1/2000, 1/2000];
%arrival_rate = 0.995*[1/2000, 1/100, 1/200, 1/2000, 1/500, 1/2000, 1/20, 1/1000, 1/200, 1/2000];

