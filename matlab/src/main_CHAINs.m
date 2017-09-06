%% Main Function For Throughput-Optimal CHAIN
clear;

%% Initialization
config.network_5nodes_1Chain;

currentT = 0;
oldT = 0;

WIFI_nodes = cell(N, 1);
next_arrival_time = zeros(N, 1);
state_vec = ones(N, 1);



%% Parameters
Npoints = 3000000;
contention_time = 0; 
q_threshold = 1; % threshold for reactivation
qn = zeros(N,1);
qn_history = zeros(N, Npoints);
frame_timestamps = zeros(1, Npoints);

for i = 1:N_CHAIN_node
    WIFI_nodes{i} = wifi_node(1,1); 
end
for i = (N_CHAIN_node+1):N
    WIFI_nodes{i} = wifi_node(0,1);
end

%% Parameters
MAX_TOL = 1e-10;


%% Main Program
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
            if winner_id <= N_CHAIN_node
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
            if winner_id <= N_CHAIN_node
                % The winner supports CHAIN
                for m=1:N_CHAIN_group
                    if ismember(winner_id, CHAINs{m})
                        for j=1:length(CHAINs{m})
                            if qn(CHAINs{m}(j)) >= q_threshold
                                qn(CHAINs{m}(j)) = qn(CHAINs{m}(j)) - 1;  % service is deterministic
                                total_TX_time = total_TX_time + Tpkt;
                                round_idx = round_idx + power(2,CHAINs{m}(j)-1);
                            end
                        end
                        break
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


