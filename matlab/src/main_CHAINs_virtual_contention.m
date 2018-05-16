%% Main Function For Throughput-Optimal CHAIN
clear;
clc;
tic;
%% Initialization
%config_virtual_contention.network_6nodes_1Chain_1DCF;
%config_virtual_contention.network_4nodes_1Chain_1DCF;
%config_virtual_contention.network_3nodes_1Chain_1DCF;
%config_virtual_contention.network_3nodes_1Chain_0DCF;
%config_virtual_contention.network_5nodes_1Chain_2DCF;
%config_virtual_contention.network_5nodes_1Chain_0DCF;
%config_virtual_contention.network_10nodes_1Chain_0DCF;
%config_virtual_contention.network_20nodes_1Chain_0DCF_2groups;
config_virtual_contention.network_30nodes_1Chain_0DCF_3groups;

currentT = 0;
oldT = 0;

WIFI_nodes = cell(N, 1);
next_arrival_time = zeros(N, 1);
state_vec = zeros(N, 1);
CHAIN_head_id = 0;
CHAIN_tail_id = 0;
CHAIN_low_head_id = 0;
CHAIN_low_tail_id = 0;
CHAIN_high_head_id = 0;
CHAIN_high_tail_id = 0;
Tput_shortterm = zeros(N, 1); % monitor the short-term empirical throughput
time_since_last_CHAIN_update = 0;

total_contention_interval = 0;
total_contention_time = 0;
queue_length_calculus = zeros(N,1);
time_since_last_queue_update = 0;
dummy_packet_count = zeros(N,1);

%% Parameters
Npoints = 3000000;
contention_time = 0; 
qn = zeros(N,1);
qn_history = zeros(N, Npoints);
frame_timestamps = zeros(1, Npoints);

for i = 1:N_CHAIN_node
    WIFI_nodes{i} = wifi_node(1,0,CW_Min,CW_Max); 
end
for i = (N_CHAIN_node+1):N
    WIFI_nodes{i} = wifi_node(0,1,CW_Min,CW_Max);
    state_vec(i) = 1;
end

%% Parameters
MAX_TOL = 1e-10;


%% Main Program
count = 0;
round_count = zeros(2^N, 1);
for i=1:N
    next_arrival_time(i) = get_interarrival_time(arrival_rate(i),number_of_mini_slots_per_unit_time);
end
while currentT < simT
    count = count + 1;
    total_TX_time = 0;
    round_idx = 1;

%% Original CHAIN
% Chain may break due to a link with an empty queue
% For simplicity, we consider that there is only one CHAIN (by choosing a small delta in the CHAIN paper)
% For fair comparison, we shall consider that CHAIN is closed
    if strcmp(Mode,'Original-CHAIN')
        % Contention
        winner_id = -1;
        if (currentT - time_since_last_CHAIN_update) >= Tupdate || currentT == 0
           % Update CHAIN piggyback relations
           % TODO!!!
           [val_sorted, id_sorted] = sort(Tput_shortterm(1:N_CHAIN_node), 'descend');
           CHAIN_head_id = id_sorted(1);
           CHAIN_tail_id = id_sorted(N_CHAIN_node);
           for m=1:N_CHAIN_node
               if m == 1 
                   WIFI_nodes{id_sorted(m)}.succ_id = id_sorted(m+1);
                   WIFI_nodes{id_sorted(m)}.pred_id = id_sorted(N_CHAIN_node);
               else
                   if m == N_CHAIN_node
                       WIFI_nodes{id_sorted(m)}.succ_id = id_sorted(1);
                       WIFI_nodes{id_sorted(m)}.pred_id = id_sorted(m-1);
                   else
                       WIFI_nodes{id_sorted(m)}.succ_id = id_sorted(m+1);
                       WIFI_nodes{id_sorted(m)}.pred_id = id_sorted(m-1);                      
                   end
               end
           end
           Tput_shortterm = zeros(N,1);
           time_since_last_CHAIN_update = currentT;
        end
        contention_time = get_contention_time(Tcont);
        % Each node is active if it has at least one packetis active if
        contention_set = find(qn >= 1);
        if ~isempty(contention_set)
            total_contention_interval = total_contention_interval + 1;
            winner_id = contention_set(randi([1, length(contention_set)]));
            currentT = currentT + contention_time;
            total_contention_time = total_contention_time + contention_time;
            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
            time_since_last_queue_update = currentT;
        end    
        % Services
        if winner_id > 0
            if winner_id <= N_CHAIN_node
                % The winner supports CHAIN
                currentT = currentT + Tpkt; 
                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                time_since_last_queue_update = currentT;
                qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                Tput_shortterm(winner_id) = Tput_shortterm(winner_id) + 1;
                round_idx = round_idx + power(2,winner_id-1); 
                next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                while next_in_CHAIN ~= winner_id
                    if qn(next_in_CHAIN) > 0
                        currentT = currentT + Tpkt; 
                        queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                        time_since_last_queue_update = currentT;
                        qn(next_in_CHAIN) = qn(next_in_CHAIN) - 1;  % service is deterministic
                        Tput_shortterm(next_in_CHAIN) = Tput_shortterm(next_in_CHAIN) + 1;
                        round_idx = round_idx + power(2,next_in_CHAIN-1);
                        next_in_CHAIN = WIFI_nodes{next_in_CHAIN}.succ_id;
                    else
                        % CHAIN breaks and the next contention interval
                        % starts
                        break
                    end
                end
            else
                % The winner follows DCF
                currentT = currentT + Tpkt; 
                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                time_since_last_queue_update = currentT;
                qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                round_idx = round_idx + power(2,winner_id-1);                
            end
        else
            if winner_id == -1
                % there is no packet available in the network
                currentT = min(next_arrival_time);
                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                time_since_last_queue_update = currentT;
            end            
        end
    end
    
%% Original CHAIN with two chains
% Chain may break due to a link with an empty queue
% For simplicity, we consider that there are possibly two CHAINs (by choosing a proper delta in the CHAIN paper)
% For fair comparison, we shall consider that CHAIN is closed
    if strcmp(Mode,'Original-CHAIN-two-chains')
        % Contention
        winner_id = -1;
        if (currentT - time_since_last_CHAIN_update) >= Tupdate || currentT == 0
           % Update CHAIN piggyback relations
           % TODO!!!
           [val_sorted, id_sorted] = sort(Tput_shortterm(1:N_CHAIN_node), 'descend');           
           id_below_delta = id_sorted(find(Tput_shortterm(id_sorted) < time_since_last_CHAIN_update*Tput_shortterm_delta));
           id_above_delta = setdiff(id_sorted, id_below_delta, 'stable'); % 'stable' mean no sorting after setdiff
           if ~isempty(id_below_delta)
               CHAIN_low_head_id = id_below_delta(1);
               CHAIN_low_tail_id = id_below_delta(length(id_below_delta));
           end
           if ~isempty(id_above_delta)
               CHAIN_high_head_id = id_above_delta(1);
               CHAIN_high_tail_id = id_above_delta(length(id_above_delta)); 
           end
           for m=1:length(id_above_delta)
               if mod(m+1,length(id_above_delta)) ~= 0
                   WIFI_nodes{id_above_delta(m)}.succ_id = id_above_delta(mod(m+1,length(id_above_delta)));
               else
                   WIFI_nodes{id_above_delta(m)}.succ_id = id_above_delta(length(id_above_delta));
               end
               if mod(m-1,length(id_above_delta)) ~= 0
                   WIFI_nodes{id_above_delta(m)}.pred_id = id_above_delta(mod(m-1,length(id_above_delta))); 
               else
                   WIFI_nodes{id_above_delta(m)}.pred_id = id_above_delta(length(id_above_delta));
               end
           end
           for m=1:length(id_below_delta)
               if mod(m+1,length(id_below_delta)) ~= 0
                   WIFI_nodes{id_below_delta(m)}.succ_id = id_below_delta(mod(m+1,length(id_below_delta)));
               else
                   WIFI_nodes{id_below_delta(m)}.succ_id = id_below_delta(length(id_below_delta));
               end
               if mod(m-1,length(id_below_delta)) ~= 0
                   WIFI_nodes{id_below_delta(m)}.pred_id = id_below_delta(mod(m-1,length(id_below_delta))); 
               else
                   WIFI_nodes{id_below_delta(m)}.pred_id = id_below_delta(length(id_below_delta));
               end
           end           
           Tput_shortterm = zeros(N,1);
           time_since_last_CHAIN_update = currentT;
        end
        contention_time = get_contention_time(Tcont);
        % Each node is active if it has at least one packetis active if
        contention_set = find(qn >= 1);
        if ~isempty(contention_set)
            total_contention_interval = total_contention_interval + 1;
            currentT = currentT + contention_time;
            total_contention_time = total_contention_time + contention_time;
            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
            time_since_last_queue_update = currentT;
            winner_id = contention_set(randi([1, length(contention_set)]));
        end    
        % Services
        if winner_id > 0
            if winner_id <= N_CHAIN_node
                % The winner supports CHAIN
                currentT = currentT + Tpkt; 
                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                time_since_last_queue_update = currentT;
                qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                Tput_shortterm(winner_id) = Tput_shortterm(winner_id) + 1;
                round_idx = round_idx + power(2,winner_id-1); 
                next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                while next_in_CHAIN ~= winner_id
                    if qn(next_in_CHAIN) > 0
                        currentT = currentT + Tpkt; 
                        queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                        time_since_last_queue_update = currentT;
                        qn(next_in_CHAIN) = qn(next_in_CHAIN) - 1;  % service is deterministic
                        Tput_shortterm(next_in_CHAIN) = Tput_shortterm(next_in_CHAIN) + 1;
                        round_idx = round_idx + power(2,next_in_CHAIN-1);
                        next_in_CHAIN = WIFI_nodes{next_in_CHAIN}.succ_id;
                    else
                        % CHAIN breaks and the next contention interval
                        % starts
                        break
                    end
                end
            else
                % The winner follows DCF
                currentT = currentT + Tpkt; 
                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                time_since_last_queue_update = currentT;
                qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                round_idx = round_idx + power(2,winner_id-1);                
            end
        else
            if winner_id == -1
                % there is no packet available in the network
                currentT = min(next_arrival_time);
                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                time_since_last_queue_update = currentT;
            end
        end
    end
    
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
            total_contention_interval = total_contention_interval + 1;
            currentT = currentT + contention_time;
            total_contention_time = total_contention_time + contention_time;
            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
            time_since_last_queue_update = currentT;
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
                    currentT = currentT + Tpkt; 
                    queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                    time_since_last_queue_update = currentT;                    
                    qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                    round_idx = round_idx + power(2,winner_id-1);
                    next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                    while next_in_CHAIN ~= winner_id
                        if qn(next_in_CHAIN) > 0
                            currentT = currentT + Tpkt; 
                            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                            time_since_last_queue_update = currentT;
                            qn(next_in_CHAIN) = qn(next_in_CHAIN) - 1;  % service is deterministic
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
                    currentT = currentT + Tpkt; 
                    queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                    time_since_last_queue_update = currentT;
                    qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                    round_idx = round_idx + power(2,winner_id-1);
                    
                    % The rest of CHAIN can piggyback with possibly dummy
                    % packets
                    next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                    while next_in_CHAIN ~= winner_id
                        if qn(next_in_CHAIN) > 0
                            currentT = currentT + Tpkt; 
                            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                            time_since_last_queue_update = currentT;
                            qn(next_in_CHAIN) = qn(next_in_CHAIN) - 1;  % service is deterministic
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        else
                            % CHAIN shall not break and can transmit dummy
                            % packets
                            currentT = currentT + Tdummy; 
                            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                            time_since_last_queue_update = currentT;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        end
                        next_in_CHAIN = WIFI_nodes{next_in_CHAIN}.succ_id;
                    end                    
                end
            else
                % The winner follows DCF
                currentT = currentT + Tpkt; 
                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                time_since_last_queue_update = currentT;
                qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
                round_idx = round_idx + power(2,winner_id-1);
            end  
        else
            if winner_id == -1
                % there is no packet available in the network
                currentT = min(next_arrival_time);
                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                time_since_last_queue_update = currentT;
            end
        end
    end

%{
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
%}        

%{    
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
%}

%{    
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
%}
   
%{    
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
%}

%{    
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
%}
    
%% Conventional Pure-DCF mode
    if strcmp(Mode,'Pure-DCF')
        % Contention
        winner_id = -1;
        contention_time = get_contention_time(Tcont);
        % Each node is active if it has at least one packet
        contention_set = find(qn >= 1);
        if ~isempty(contention_set)
            total_contention_interval = total_contention_interval + 1;
            currentT = currentT + contention_time;
            total_contention_time = total_contention_time + contention_time;
            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
            time_since_last_queue_update = currentT;
            winner_id = contention_set(randi([1, length(contention_set)]));
        end 
        if winner_id > 0
            % Time evolution
            currentT = currentT + Tpkt; 
            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
            time_since_last_queue_update = currentT;
            
        	qn(winner_id) = qn(winner_id) - 1;  % service is deterministic
            round_idx = round_idx + power(2,winner_id-1);
        else
            if winner_id == -1
                % there is no packet available in the network
                currentT = min(next_arrival_time);
                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                time_since_last_queue_update = currentT;
            end           
        end
    end
    
%% Common parts    
    % Arrivals
    for i=1:N
        while next_arrival_time(i) <= currentT + MAX_TOL
            qn(i) = qn(i) + 1;
            next_arrival_time(i) = next_arrival_time(i) + get_interarrival_time(arrival_rate(i),number_of_mini_slots_per_unit_time);
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
average_contention_time_per_interval = total_contention_time/total_contention_interval;

%% Plotting
%plot(frame_timestamps(1:count), qn_history(:,1:count));
createfigure(frame_timestamps(1:count), qn_history(:,1:count));
average_qn = queue_length_calculus/currentT;
for i=1:N
fprintf('average queue length: %.2f\n', average_qn(i));
end
fprintf('total average queue length: %.2f\n', sum(average_qn));
fprintf('total average delay: %.2f (ms)\n', sum(average_qn)/sum(arrival_rate));
fprintf('average contention time per interval: %.3f (ms)\n', average_contention_time_per_interval);
fprintf('contention time percentage: %.3f %%\n', total_contention_time/currentT*100);
toc;

