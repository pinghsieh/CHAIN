%% Main Function For Throughput-Optimal CHAIN
clear;
clc;
tic;

%% Part 1. Initialization
simT = 1000000;
currentT = 0;
oldT = 0;
Tpkt = 1;
Tyield = 1;
Tcont = 1;

Mode = 'Qth-based';
%Mode = 'Qth-cross-piggyback';
%Mode = 'Qth-plus-Contention';
N_CHAIN = 3;
N_DCF = 0;
N = N_CHAIN + N_DCF;
is_CHAIN = [ones(N_CHAIN,1); zeros(N_DCF,1)];
is_DCF = [zeros(N_CHAIN,1); ones(N_DCF,1)];
round_count = zeros(2^N, 1);

WIFI_nodes = cell(N, 1);
next_arrival_time = zeros(N, 1);
state_vec = ones(N, 1);


% Arrival rates for Qth-based
%arrival_rate = 0.96*[1/15*0.3, 1/10*0.2, 1/15*0.4, 1/15*0.05, 1/15*0.05];
%arrival_rate = 0.96*[1/15*0.15, 1/10*0.1, 1/15*0.4, 1/15*0.2, 1/15*0.15];
%arrival_rate = 1*[0.04, 0.01, 0.33, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01];
arrival_rate = 1*[0.04, 0.02, 0.2, 0.2, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01];
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


% Arrival rates for Qth-cross-piggyback



% Parameters
Npoints = 3000000;

contention_time = 0; 
q_threshold = 10; % threshold for reactivation
qn = zeros(N,1);
qn_history = zeros(N, Npoints);
frame_timestamps = zeros(1, Npoints);

for i = 1:N_CHAIN
    WIFI_nodes{i} = wifi_node(1,1); 
end
for i = (N_CHAIN+1):N
    WIFI_nodes{i} = wifi_node(0,1);
end

%% Parameters
MAX_TOL = 1e-10;


%% Part 2: Main Program
count = 0;
for i=1:N
    next_arrival_time(i) = get_interarrival_time(arrival_rate(i));
end
while currentT < simT
    count = count + 1;
    total_TX_time = 0;
    round_idx = 1;

% Q-threshold plus contention-based re-admission
    if strcmp(Mode,'Qth-plus-Contention')
        % Contention
        winner_id = -1;
        contention_time = get_contention_time(Tcont);
        inactive_set = find(state_vec == 0);
        active_set = find(state_vec == 1);
        contention_set_part1 = [];
        contention_set_part2 = [];
        if ~isempty(inactive_set)
            contention_set_part1 = find(qn.*(state_vec == 0) >= q_threshold);
        end
        if ~isempty(active_set)
            contention_set_part2 = find(qn.*(state_vec == 1) > 0);
        end
        contention_set = [contention_set_part1; contention_set_part2];
    
        if ~isempty(contention_set)
            winner_id = contention_set(randi([1, length(contention_set)]));
        end
    
        % Services
        if winner_id > 0
            if winner_id <= N_CHAIN
                % The winner supports CHAIN
                if WIFI_nodes{winner_id}.is_active == 1
                    for m=0:N_CHAIN-1
                        j = mod(winner_id+m-1,N_CHAIN) + 1;
                        if WIFI_nodes{j}.is_active == 1
                            if qn(j) > 0
                                qn(j) = qn(j) - 1;  % service is deterministic
                                total_TX_time = total_TX_time + Tpkt;
                                round_idx = round_idx + power(2,j-1);
                            else
                                WIFI_nodes{j}.is_active = 0;
                                state_vec(j) = 0;
                                %total_TX_time = total_TX_time + Tyield;
                                break
                            end
                        end
                    end
                else
                    WIFI_nodes{winner_id}.is_active = 1;
                    state_vec(winner_id) = 1;
                    qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                    total_TX_time = total_TX_time + Tpkt;
                    round_idx = round_idx + power(2,winner_id-1);
                end
            else
                % The winner follows DCF
                qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                total_TX_time = total_TX_time + Tpkt;
                round_idx = round_idx + power(2,winner_id-1);
            end
            % Time evolution
            total_frame_time = contention_time + total_TX_time;
            oldT = currentT;
            currentT = currentT + total_frame_time;  
        else
            oldT = currentT;
            currentT = min(next_arrival_time);
        end
    end
    
% Q-threshold based  
    if strcmp(Mode,'Qth-based')
        % Contention
        winner_id = -1;
        contention_time = get_contention_time(Tcont);
        contention_set_part1 = find(qn.*is_CHAIN >= q_threshold);
        contention_set_part2 = find(qn.*is_DCF > 0);
        contention_set = [contention_set_part1; contention_set_part2];
    
        if ~isempty(contention_set)
            winner_id = contention_set(randi([1, length(contention_set)]));
        end
        % Services
        if winner_id > 0
            if winner_id <= N_CHAIN
                % The winner supports CHAIN
                for j=1:N_CHAIN
                    if qn(j) >= q_threshold
                        qn(j) = qn(j) - 1;  % service is deterministic
                        total_TX_time = total_TX_time + Tpkt;
                        round_idx = round_idx + power(2,j-1);
                    end
                end
            else
                % The winner follows DCF
                qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                total_TX_time = total_TX_time + Tpkt;
                round_idx = round_idx + power(2,winner_id-1);
            end
            % Time evolution
            total_frame_time = contention_time + total_TX_time;
            oldT = currentT;
            currentT = currentT + total_frame_time;  
        else
            oldT = currentT;
            currentT = min(next_arrival_time);
        end
    end              

% Q-threshold with cross piggybacking  
    if strcmp(Mode,'Qth-cross-piggyback')
        % Contention
        winner_id = -1;
        contention_time = get_contention_time(Tcont);
        contention_set_part1 = find(qn.*is_CHAIN >= q_threshold);
        contention_set_part2 = find(qn.*is_DCF > 0);
        contention_set = [contention_set_part1; contention_set_part2];
    
        if ~isempty(contention_set)
            winner_id = contention_set(randi([1, length(contention_set)]));
        end
        % Services
        if winner_id > 0
            if winner_id <= N_CHAIN
                % The winner supports CHAIN
                for j=1:N_CHAIN
                    if qn(j) >= q_threshold
                        qn(j) = qn(j) - 1;  % service is deterministic
                        total_TX_time = total_TX_time + Tpkt;
                        round_idx = round_idx + power(2,j-1);
                    end
                end
            else
                % The winner follows DCF
                qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                total_TX_time = total_TX_time + Tpkt;
                round_idx = round_idx + power(2,winner_id-1);
                % CHAIN nodes can piggyback a DCF transmission
                for j=1:N_CHAIN
                    if qn(j) >= q_threshold
                        qn(j) = qn(j) - 1;  % service is deterministic
                        total_TX_time = total_TX_time + Tpkt;
                        round_idx = round_idx + power(2,j-1);
                    end
                end
            end
            % Time evolution
            total_frame_time = contention_time + total_TX_time;
            oldT = currentT;
            currentT = currentT + total_frame_time;  
        else
            oldT = currentT;
            currentT = min(next_arrival_time);
        end
    end 
    
%% Common parts    
    % Arrivals
    for i=1:N
        while next_arrival_time(i) <= currentT + MAX_TOL
            qn(i) = qn(i) + 1;
            next_arrival_time(i) = next_arrival_time(i) + get_interarrival_time(arrival_rate(i));
        end
    end
    
    % Round Counting
    round_count(round_idx,1) = round_count(round_idx,1) + 1;
    
    % Update history
    qn_history(:,count) = qn;
    frame_timestamps(count) = currentT;
end

%% Post-processing
round_count_normalize = round_count/(sum(round_count));

%% Plotting
%plot(frame_timestamps(1:count), qn_history(:,1:count));
createfigure(frame_timestamps(1:count), qn_history(:,1:count));

toc;
