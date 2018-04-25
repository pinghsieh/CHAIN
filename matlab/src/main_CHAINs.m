%% Main Function For Throughput-Optimal CHAIN
clear;
clc;
tic;
%% Initialization
%config.network_5nodes_1Chain;
%config.network_6nodes_1Chain_1DCF;
config.network_4nodes_1Chain_1DCF;
%config.network_3nodes_1Chain_1DCF;
%config.network_3nodes_1Chain_0DCF;
%config.network_5nodes_1Chain_2DCF;

currentT = 0;
oldT = 0;

WIFI_nodes = cell(N, 1);
next_arrival_time = zeros(N, 1);
state_vec = zeros(N, 1);
CHAIN_head_id = 0;
CHAIN_tail_id = 0;


%% Parameters
Npoints = 3000000;
contention_time = 0; 
qn = zeros(N,1);
qn_history = zeros(N, Npoints);
frame_timestamps = zeros(1, Npoints);

for i = 1:N_CHAIN_node
    WIFI_nodes{i} = wifi_node(1,0); 
end
for i = (N_CHAIN_node+1):N
    WIFI_nodes{i} = wifi_node(0,1);
    state_vec(i) = 1;
end

%% Parameters
MAX_TOL = 1e-10;


%% Main Program
count = 0;
round_count = zeros(2^N, 1);
for i=1:N
    next_arrival_time(i) = get_interarrival_time(arrival_rate(i));
end
while currentT < simT
    count = count + 1;
    total_TX_time = 0;
    round_idx = 1;

%% Q-threshold plus contention-based re-admission
% Chain may break due to a link with an empty queue
    if strcmp(Mode,'Qth-plus-Contention')
        % Contention
        winner_id = -1;
        contention_time = get_contention_time(Tcont);
        % CHAIN node is active if
        % (1) it is already in CHAIN and has at least one packet
        % (2) it it currently not in CHAIN but has more than or equal to
        % q_threshold packets
        contention_set_part1 = find(qn(1:N_CHAIN_node).*(state_vec(1:N_CHAIN_node) == 0) >= q_threshold);
        contention_set_part2 = find(qn(1:N_CHAIN_node).*(state_vec(1:N_CHAIN_node) == 1) >= 1);
        % DCF node is active if it has at least one packet
        contention_set_part3 = find(qn.*is_DCF >= 1);
        contention_set = [contention_set_part1; contention_set_part2; contention_set_part3];
        if ~isempty(contention_set)
            winner_id = contention_set(randi([1, length(contention_set)]));
        end
    
        % Services
        if winner_id > 0
            if winner_id <= N_CHAIN_node
                % The winner supports CHAIN
                if WIFI_nodes{winner_id}.is_active == 1
                    % The winner is already in CHAIN
                    % TODO: In this case, CHAIN might break if a link in CHAIN
                    % has an empty queue
                    qn(winner_id) = max(0, qn(winner_id) - 1);  % service is deterministic
                    total_TX_time = total_TX_time + Tpkt;
                    round_idx = round_idx + power(2,winner_id-1);
                    next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                    while next_in_CHAIN ~= winner_id
                        if qn(next_in_CHAIN) > 0
                            qn(next_in_CHAIN) = qn(next_in_CHAIN) - 1;  % service is deterministic
                            total_TX_time = total_TX_time + Tpkt;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        else
                            % CHAIN breaks and needs to update predecessors
                            % and successors
                            WIFI_nodes{next_in_CHAIN}.is_active = 0;
                            state_vec(next_in_CHAIN) = 0;
                            % Update CHAIN head or tail
                            if next_in_CHAIN == CHAIN_tail_id
                                CHAIN_tail_id = WIFI_nodes{next_in_CHAIN}.pred_id;
                            else
                                if next_in_CHAIN == CHAIN_head_id
                                    CHAIN_head_id = WIFI_nodes{next_in_CHAIN}.succ_id;       
                                end
                            end
                            WIFI_nodes{WIFI_nodes{next_in_CHAIN}.succ_id}.pred_id = WIFI_nodes{next_in_CHAIN}.pred_id;
                            WIFI_nodes{WIFI_nodes{next_in_CHAIN}.pred_id}.succ_id = WIFI_nodes{next_in_CHAIN}.succ_id;
                            WIFI_nodes{next_in_CHAIN}.succ_id = 0;
                            WIFI_nodes{next_in_CHAIN}.pred_id = 0;
                            break
                        end
                        next_in_CHAIN = WIFI_nodes{next_in_CHAIN}.succ_id;
                    end
                else
                    % The winner is currently not in CHAIN
                    % TODO: The winner joins the CHAIN as the head, and the
                    % rest of CHAIN can piggyback
                    WIFI_nodes{winner_id}.is_active = 1;
                    if CHAIN_tail_id ~= 0
                        WIFI_nodes{winner_id}.pred_id = CHAIN_tail_id;
                        WIFI_nodes{CHAIN_tail_id}.succ_id = winner_id;
                    else
                        % if there was no link in CHAIN
                        WIFI_nodes{winner_id}.pred_id = winner_id;
                        CHAIN_tail_id = winner_id;
                    end
                    if CHAIN_head_id ~= 0
                        WIFI_nodes{winner_id}.succ_id = CHAIN_head_id;
                        WIFI_nodes{CHAIN_head_id}.pred_id = winner_id;
                    else
                        % if there was no link in CHAIN
                        WIFI_nodes{winner_id}.succ_id = winner_id;
                    end
                    CHAIN_head_id = winner_id;                  
                    state_vec(winner_id) = 1;
                    qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                    total_TX_time = total_TX_time + Tpkt;
                    round_idx = round_idx + power(2,winner_id-1);
                    
                    % The rest of CHAIN can piggyback with possibly dummy
                    % packets
                    next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                    while next_in_CHAIN ~= winner_id
                        if qn(next_in_CHAIN) > 0
                            qn(next_in_CHAIN) = qn(next_in_CHAIN) - 1;  % service is deterministic
                            total_TX_time = total_TX_time + Tpkt;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        else
                            % CHAIN shall not break and can transmit dummy
                            % packets
                            total_TX_time = total_TX_time + Tdummy;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        end
                        next_in_CHAIN = WIFI_nodes{next_in_CHAIN}.succ_id;
                    end                    
                end
            else
                % The winner follows DCF
                qn(winner_id) = max(0, qn(winner_id) - 1);  % service is deterministic
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

%% Q-threshold plus contention-based re-admission, and DCF nodes are given higher priority
% Chain may break due to a link with an empty queue
    if strcmp(Mode,'Qth-plus-Contention-DCF-priority')
        % Contention
        winner_id = -1;
        contention_time = get_contention_time(Tcont);
        % CHAIN node is active if
        % (1) it is already in CHAIN and has at least one packet
        % (2) it it currently not in CHAIN but has more than or equal to
        % q_threshold packets
        contention_set_part1 = find(qn(1:N_CHAIN_node).*(state_vec(1:N_CHAIN_node) == 0) >= q_threshold);
        contention_set_part2 = find(qn(1:N_CHAIN_node).*(state_vec(1:N_CHAIN_node) == 1) >= 1);
        % DCF node is active if it has at least one packet
        contention_set_part3 = find(qn.*is_DCF >= 1);
        contention_set = [contention_set_part1; contention_set_part3];
        if ~isempty(contention_set)
            winner_id = contention_set(randi([1, length(contention_set)]));
        else
            if ~isempty(contention_set_part2)
                 winner_id = contention_set_part2(randi([1, length(contention_set_part2)]));
            end
        end
    
        % Services
        if winner_id > 0
            if winner_id <= N_CHAIN_node
                % The winner supports CHAIN
                if WIFI_nodes{winner_id}.is_active == 1
                    % The winner is already in CHAIN
                    % TODO: In this case, CHAIN might break if a link in CHAIN
                    % has an empty queue
                    qn(winner_id) = max(0, qn(winner_id) - 1);  % service is deterministic
                    total_TX_time = total_TX_time + Tpkt;
                    round_idx = round_idx + power(2,winner_id-1);
                    next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                    while next_in_CHAIN ~= winner_id
                        if qn(next_in_CHAIN) > 0
                            qn(next_in_CHAIN) = qn(next_in_CHAIN) - 1;  % service is deterministic
                            total_TX_time = total_TX_time + Tpkt;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        else
                            % CHAIN breaks and needs to update predecessors
                            % and successors
                            WIFI_nodes{next_in_CHAIN}.is_active = 0;
                            state_vec(next_in_CHAIN) = 0;
                            % Update CHAIN head or tail
                            if next_in_CHAIN == CHAIN_tail_id
                                CHAIN_tail_id = WIFI_nodes{next_in_CHAIN}.pred_id;
                            else
                                if next_in_CHAIN == CHAIN_head_id
                                    CHAIN_head_id = WIFI_nodes{next_in_CHAIN}.succ_id;       
                                end
                            end
                            WIFI_nodes{WIFI_nodes{next_in_CHAIN}.succ_id}.pred_id = WIFI_nodes{next_in_CHAIN}.pred_id;
                            WIFI_nodes{WIFI_nodes{next_in_CHAIN}.pred_id}.succ_id = WIFI_nodes{next_in_CHAIN}.succ_id;
                            WIFI_nodes{next_in_CHAIN}.succ_id = 0;
                            WIFI_nodes{next_in_CHAIN}.pred_id = 0;
                            break
                        end
                        next_in_CHAIN = WIFI_nodes{next_in_CHAIN}.succ_id;
                    end
                else
                    % The winner is currently not in CHAIN
                    % TODO: The winner joins the CHAIN as the head, and the
                    % rest of CHAIN can piggyback
                    WIFI_nodes{winner_id}.is_active = 1;
                    if CHAIN_tail_id ~= 0
                        WIFI_nodes{winner_id}.pred_id = CHAIN_tail_id;
                        WIFI_nodes{CHAIN_tail_id}.succ_id = winner_id;
                    else
                        % if there was no link in CHAIN
                        WIFI_nodes{winner_id}.pred_id = winner_id;
                        CHAIN_tail_id = winner_id;
                    end
                    if CHAIN_head_id ~= 0
                        WIFI_nodes{winner_id}.succ_id = CHAIN_head_id;
                        WIFI_nodes{CHAIN_head_id}.pred_id = winner_id;
                    else
                        % if there was no link in CHAIN
                        WIFI_nodes{winner_id}.succ_id = winner_id;
                    end
                    CHAIN_head_id = winner_id;                  
                    state_vec(winner_id) = 1;
                    qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                    total_TX_time = total_TX_time + Tpkt;
                    round_idx = round_idx + power(2,winner_id-1);
                    
                    % The rest of CHAIN can piggyback with possibly dummy
                    % packets
                    next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                    while next_in_CHAIN ~= winner_id
                        if qn(next_in_CHAIN) > 0
                            qn(next_in_CHAIN) = qn(next_in_CHAIN) - 1;  % service is deterministic
                            total_TX_time = total_TX_time + Tpkt;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        else
                            % CHAIN shall not break and can transmit dummy
                            % packets
                            total_TX_time = total_TX_time + Tdummy;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        end
                        next_in_CHAIN = WIFI_nodes{next_in_CHAIN}.succ_id;
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
        

    
%% Q-threshold based 
% If an active CHAIN node wins contention, then all the active CHAIN nodes
% can transmit one packet in this interval
% A CHAIN node is active if
% (1) it is already active and has at least one packet
% (2) it it currently not active but has more than or equal to
% q_threshold packets at the beginning of the current interval

    if strcmp(Mode,'Qth-based')
        % Contention
        winner_id = -1;
        contention_time = get_contention_time(Tcont);
        % CHAIN node is active if
        % (1) it is already active and has at least one packet
        % (2) it it currently not active but has more than or equal to
        % q_threshold packets
        state_vec(1:N_CHAIN_node) = state_vec(1:N_CHAIN_node) - ((state_vec(1:N_CHAIN_node) == 1).*(qn(1:N_CHAIN_node) == 0));
        state_vec(1:N_CHAIN_node) = ((state_vec(1:N_CHAIN_node) + (qn(1:N_CHAIN_node) >= q_threshold)) > 0.5);
        
        contention_set_part1 = find(state_vec(1:N_CHAIN_node) == 1);
        % DCF node is active if it has at least one packet
        contention_set_part2 = find(qn.*is_DCF >= 1);
        contention_set = [contention_set_part1; contention_set_part2];
        if ~isempty(contention_set)
            winner_id = contention_set(randi([1, length(contention_set)]));
        end
        
        % Services
        if winner_id > 0
            if winner_id <= N_CHAIN_node
                % The winner supports CHAIN
                % TODO: In this case, an active CHAIN node with a non-empty
                % queue can transmit.
                for j=1:N_CHAIN_node
                    if state_vec(j) == 1 
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

%% Q-threshold based and let DCF nodes have higher priority
% During contention, all DCF nodes have higher priority than CHAIN nodes
% If an active CHAIN node wins contention, then all the active CHAIN nodes
% can transmit one packet in this interval
% A CHAIN node is active if
% (1) it is already active and has at least one packet
% (2) it it currently not active but has more than or equal to
% q_threshold packets at the beginning of the current interval

    if strcmp(Mode,'Qth-based-DCF-priority')
        % Contention
        winner_id = -1;
        contention_time = get_contention_time(Tcont);
        % CHAIN node is active if
        % (1) it is already active and has at least one packet
        % (2) it it currently not active but has more than or equal to
        % q_threshold packets
        state_vec(1:N_CHAIN_node) = state_vec(1:N_CHAIN_node) - ((state_vec(1:N_CHAIN_node) == 1).*(qn(1:N_CHAIN_node) == 0));
        state_vec(1:N_CHAIN_node) = ((state_vec(1:N_CHAIN_node) + (qn(1:N_CHAIN_node) >= q_threshold)) > 0.5);
        
        contention_set_part1 = find(state_vec(1:N_CHAIN_node) == 1);
        % DCF node is active if it has at least one packet
        contention_set_part2 = find(qn.*is_DCF >= 1);
        if ~isempty(contention_set_part2)
            winner_id = contention_set_part2(randi([1, length(contention_set_part2)]));
        else
            if ~isempty(contention_set_part1)
                winner_id = contention_set_part1(1);
            end
        end
        
        % Services
        if winner_id > 0
            if winner_id <= N_CHAIN_node
                % The winner supports CHAIN
                % TODO: In this case, an active CHAIN node with a non-empty
                % queue can transmit.
                for j=1:N_CHAIN_node
                    if state_vec(j) == 1 
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

%% Q-threshold based and use Chain-wide contention 
% If an active CHAIN node wins contention, then all the active CHAIN nodes
% can transmit one packet in this interval
% A CHAIN node is active if
% (1) it is already active and has at least one packet
% (2) it it currently not active but has more than or equal to
% q_threshold packets at the beginning of the current interval
% Only the representer of Chain can contend for channel access

    if strcmp(Mode,'Qth-based-Chain-wide-contention')
        % Contention
        winner_id = -1;
        contention_time = get_contention_time(Tcont);
        % CHAIN node is active if
        % (1) it is already active and has at least one packet
        % (2) it it currently not active but has more than or equal to
        % q_threshold packets
        state_vec(1:N_CHAIN_node) = state_vec(1:N_CHAIN_node) - ((state_vec(1:N_CHAIN_node) == 1).*(qn(1:N_CHAIN_node) == 0));
        state_vec(1:N_CHAIN_node) = ((state_vec(1:N_CHAIN_node) + (qn(1:N_CHAIN_node) >= q_threshold)) > 0.5);
        
        contention_set_part1 = find(state_vec(1:N_CHAIN_node) == 1);
        % DCF node is active if it has at least one packet
        contention_set_part2 = find(qn.*is_DCF >= 1);
        if ~isempty(contention_set_part1)
            contention_set = [contention_set_part1(1); contention_set_part2]; 
        else
            contention_set = contention_set_part2;
        end
        if ~isempty(contention_set)
            winner_id = contention_set(randi([1, length(contention_set)]));
        end
        
        % Services
        if winner_id > 0
            if winner_id <= N_CHAIN_node
                % The winner supports CHAIN
                % TODO: In this case, an active CHAIN node with a non-empty
                % queue can transmit.
                for j=1:N_CHAIN_node
                    if state_vec(j) == 1 
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
    
%% Q-threshold with cross piggybacking
% If a DCF node wins contention, the head of CHAIN can still piggyback
    if strcmp(Mode,'Qth-cross-piggyback')
        % Contention
        winner_id = -1;
        contention_time = get_contention_time(Tcont);
        % CHAIN node is active if
        % (1) it is already in CHAIN and has at least one packet
        % (2) it it currently not in CHAIN but has more than or equal to
        % q_threshold packets
        contention_set_part1 = find(qn(1:N_CHAIN_node).*(state_vec(1:N_CHAIN_node) == 0) >= q_threshold);
        contention_set_part2 = find(qn(1:N_CHAIN_node).*(state_vec(1:N_CHAIN_node) == 1) >= 1);
        % DCF node is active if it has at least one packet
        contention_set_part3 = find(qn.*is_DCF >= 1);
        contention_set = [contention_set_part1; contention_set_part2; contention_set_part3];
        if ~isempty(contention_set)
            winner_id = contention_set(randi([1, length(contention_set)]));
        end
        
        % Services
        if winner_id > 0
            if winner_id <= N_CHAIN_node
                % The winner supports CHAIN
                if WIFI_nodes{winner_id}.is_active == 1
                    % The winner is already in CHAIN
                    % TODO: In this case, CHAIN might break if a link in CHAIN
                    % has an empty queue
                    qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                    total_TX_time = total_TX_time + Tpkt;
                    round_idx = round_idx + power(2,winner_id-1);
                    next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                    while next_in_CHAIN ~= winner_id
                        if qn(next_in_CHAIN) > 0
                            qn(next_in_CHAIN) = qn(next_in_CHAIN) - 1;  % service is deterministic
                            total_TX_time = total_TX_time + Tpkt;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        else
                            % CHAIN breaks and needs to update predecessors
                            % and successors
                            WIFI_nodes{next_in_CHAIN}.is_active = 0;
                            state_vec(next_in_CHAIN) = 0;
                            % Update CHAIN head or tail
                            if next_in_CHAIN == CHAIN_tail_id
                                CHAIN_tail_id = WIFI_nodes{next_in_CHAIN}.pred_id;
                            else
                                if next_in_CHAIN == CHAIN_head_id
                                    CHAIN_head_id = WIFI_nodes{next_in_CHAIN}.succ_id;       
                                end
                            end
                            WIFI_nodes{WIFI_nodes{next_in_CHAIN}.succ_id}.pred_id = WIFI_nodes{next_in_CHAIN}.pred_id;
                            WIFI_nodes{WIFI_nodes{next_in_CHAIN}.pred_id}.succ_id = WIFI_nodes{next_in_CHAIN}.succ_id;
                            WIFI_nodes{next_in_CHAIN}.succ_id = 0;
                            WIFI_nodes{next_in_CHAIN}.pred_id = 0;
                            break
                        end
                        next_in_CHAIN = WIFI_nodes{next_in_CHAIN}.succ_id;
                    end
                else
                    % The winner is currently not in CHAIN
                    % TODO: The winner joins the CHAIN as the head, and the
                    % rest of CHAIN can piggyback
                    WIFI_nodes{winner_id}.is_active = 1;
                    if CHAIN_tail_id ~= 0
                        WIFI_nodes{winner_id}.pred_id = CHAIN_tail_id;
                        WIFI_nodes{CHAIN_tail_id}.succ_id = winner_id;
                    else
                        % if there was no link in CHAIN
                        WIFI_nodes{winner_id}.pred_id = winner_id;
                        CHAIN_tail_id = winner_id;
                    end
                    if CHAIN_head_id ~= 0
                        WIFI_nodes{winner_id}.succ_id = CHAIN_head_id;
                        WIFI_nodes{CHAIN_head_id}.pred_id = winner_id;
                    else
                        % if there was no link in CHAIN
                        WIFI_nodes{winner_id}.succ_id = winner_id;
                    end
                    CHAIN_head_id = winner_id;                  
                    state_vec(winner_id) = 1;
                    qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                    total_TX_time = total_TX_time + Tpkt;
                    round_idx = round_idx + power(2,winner_id-1);
                    
                    % The rest of CHAIN can piggyback with possibly dummy
                    % packets
                    next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                    while next_in_CHAIN ~= winner_id
                        if qn(next_in_CHAIN) > 0
                            qn(next_in_CHAIN) = qn(next_in_CHAIN) - 1;  % service is deterministic
                            total_TX_time = total_TX_time + Tpkt;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        else
                            % CHAIN shall not break and can transmit dummy
                            % packets
                            total_TX_time = total_TX_time + Tdummy;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        end
                        next_in_CHAIN = WIFI_nodes{next_in_CHAIN}.succ_id;
                    end                    
                end
            else
                % The winner follows DCF
                qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                total_TX_time = total_TX_time + Tpkt;
                round_idx = round_idx + power(2,winner_id-1);
                
                % The head of CHAIN can piggyback and the rest of CHAIN can
                % also follow
                if CHAIN_head_id > 0
                    if qn(CHAIN_head_id) > 0
                        qn(CHAIN_head_id) = qn(CHAIN_head_id) - 1; % service is deterministic
                        total_TX_time = total_TX_time + Tpkt;
                    else
                        total_TX_time = total_TX_time + Tdummy;
                    end              
                    next_in_CHAIN = WIFI_nodes{CHAIN_head_id}.succ_id;
                    while next_in_CHAIN ~= CHAIN_head_id
                        if qn(next_in_CHAIN) > 0
                            qn(next_in_CHAIN) = qn(next_in_CHAIN) - 1;  % service is deterministic
                            total_TX_time = total_TX_time + Tpkt;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        else
                            % CHAIN breaks and needs to update predecessors
                            % and successors
                            WIFI_nodes{next_in_CHAIN}.is_active = 0;
                            state_vec(next_in_CHAIN) = 0;
                            % Update CHAIN head or tail
                            if next_in_CHAIN == CHAIN_tail_id
                                CHAIN_tail_id = WIFI_nodes{next_in_CHAIN}.pred_id;
                            else
                                if next_in_CHAIN == CHAIN_head_id
                                    CHAIN_head_id = WIFI_nodes{next_in_CHAIN}.succ_id;       
                                end
                            end
                            WIFI_nodes{WIFI_nodes{next_in_CHAIN}.succ_id}.pred_id = WIFI_nodes{next_in_CHAIN}.pred_id;
                            WIFI_nodes{WIFI_nodes{next_in_CHAIN}.pred_id}.succ_id = WIFI_nodes{next_in_CHAIN}.succ_id;
                            WIFI_nodes{next_in_CHAIN}.succ_id = 0;
                            WIFI_nodes{next_in_CHAIN}.pred_id = 0;
                            break
                        end
                        next_in_CHAIN = WIFI_nodes{next_in_CHAIN}.succ_id;
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
    
%% Conventional Pure-DCF mode
    if strcmp(Mode,'Pure-DCF')
        % Contention
        winner_id = -1;
        contention_time = get_contention_time(Tcont);
        % Each node is active if it has at least one packet
        contention_set = find(qn >= 1);
        if ~isempty(contention_set)
            winner_id = contention_set(randi([1, length(contention_set)]));
        end 
        if winner_id > 0
        	qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
            total_TX_time = total_TX_time + Tpkt;
            round_idx = round_idx + power(2,winner_id-1);
            
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
average_qn = sum(qn_history(:,1:count), 2)/count;
for i=1:N
fprintf('average queue length: %.2f\n', average_qn(i));
end
toc;

