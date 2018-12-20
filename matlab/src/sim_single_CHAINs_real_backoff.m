function [] = sim_single_CHAINs_real_backoff(network_config_code, Mode, fileID, varargin)
%% Main Function For Throughput-Optimal CHAIN
% This function requires 3 inputs, 2 optionals
%% Choose algorithm
% Qth-based: 
% Qth-cross-piggyback:
% Qth-plus-Contention:
% Qth-plus-Contention-with-dummy
% Original-CHAIN:
% Original-CHAIN-two-chains:
% Pure-DCF:

tic;
%% Check number of input arguments
if nargin >= 6
    error('myfuns:sim_single_CHAINs_real_backoff:TooManyInputs', ...
        'requires at most 3 optional inputs');
end

%% Fill in unset optional values.
switch nargin
    case 3
        rho = 0;
    case 4
        rho = varargin{1};
    case 5
        rho1 = varargin{1};
        rho2 = varargin{2};
end
%% Initialization
switch network_config_code
    case '2-1-0-0-1group'
        config_real_backoff.network_2nodes_1Chain_0DCF_0DCFsat_1group;    
    case '3-1-0-0-1group'
        config_real_backoff.network_3nodes_1Chain_0DCF_0DCFsat_1group;  
    case '10-1-0-0-1group'
        config_real_backoff.network_10nodes_1Chain_0DCF_0DCFsat_1group;   
    case '10-1-0-0-10groups-het'
        config_real_backoff.network_10nodes_1Chain_0DCF_0DCFsat_10groups_het;         
    case '500-1-0-0-1group'
        config_real_backoff.network_500nodes_1Chain_0DCF_0DCFsat_1group; 
    case '500-1-0-0-2groups'
        config_real_backoff.network_500nodes_1Chain_0DCF_0DCFsat_2groups; 
    case '500-1-0-0-2groups-het'
        config_real_backoff.network_500nodes_1Chain_0DCF_0DCFsat_2groups_het;    
    case '500-1-0-0-10groups-het'
        config_real_backoff.network_500nodes_1Chain_0DCF_0DCFsat_10groups_het;  
    case '500-1-0-0-10groups-het-unreliable'
        config_real_backoff.network_500nodes_1Chain_0DCF_0DCFsat_10groups_het_unreliable;      
    case '5-1-0'
        config_real_backoff.network_5nodes_1Chain_0DCF;   
    case '10-1-0'
        config_real_backoff.network_10nodes_1Chain_0DCF;
    case '20-1-0'
        config_real_backoff.network_20nodes_1Chain_0DCF;
    case '20-1-0-2groups'
        config_real_backoff.network_20nodes_1Chain_0DCF_2groups;
    case '20-1-0-2groups-het'
        config_real_backoff.network_20nodes_1Chain_0DCF_2groups_hetero;   
    case '20-1-10-2groups-het'
        config_real_backoff.network_20nodes_1Chain_10DCF_2groups_hetero;        
    case '20-1-0-4groups'
        config_real_backoff.network_20nodes_1Chain_0DCF_4groups; 
    case '30-1-0-10-2groups-sat'
        config_real_backoff.network_30nodes_1Chain_0DCF_10DCFsat_2groups;          
    case '30-1-0-3groups'
        config_real_backoff.network_30nodes_1Chain_0DCF_3groups; 
    case '30-1-10-3groups'
        config_real_backoff.network_30nodes_1Chain_10DCF_3groups; 
    otherwise
        config_real_backoff.network_2nodes_1Chain_0DCF_0DCFsat_1group; 
end

total_arrival = sum(arrival_rate);

%% Print log headers
switch nargin
    case 3
        fprintf(fileID, '***** Mode = %s, simT = %d (ms), rho = %.4f, total arrival = %.4f, p = %.3f, Qth = %d ***** \n', Mode, simT, rho, total_arrival, p_star, q_threshold);
    case 4
        fprintf(fileID, '***** Mode = %s, simT = %d (ms), rho = %.4f, total arrival = %.4f, p = %.3f, Qth = %d ***** \n', Mode, simT, rho, total_arrival, p_star, q_threshold);
    case 5
        fprintf(fileID, '***** Mode = %s, simT = %d (ms), rho1 = %.4f, rho2 = %.4f, total arrival = %.4f, p = %.3f, Qth = %d ***** \n', Mode, simT, rho1, rho2, total_arrival, p_star, q_threshold);
end


currentT = 0;

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
backoff = zeros(N, 1);
total_collision = 0;
total_contention_interval = 0;
total_contention_time = 0;
queue_length_calculus = zeros(N,1);
time_since_last_queue_update = 0;
dummy_packet_count = zeros(N,1);
time_last_queue_update_QLBCSMA = 0;

%% Parameters
Npoints = 500000;
Nmorepoints = 100000;
qn = zeros(N,1);
qn_slow = zeros(N,1); % for QLB-CSMA
qn_history = zeros(N, Npoints);
total_delivery_history = zeros(N, Npoints);
frame_timestamps = zeros(1, Npoints);

for i = 1:N_CHAIN_node
    WIFI_nodes{i} = wifi_node(1, 0, CW_Min, CW_Max); 
end
for i = (N_CHAIN_node+1):N
    WIFI_nodes{i} = wifi_node(0, 1, CW_Min, CW_Max);
    state_vec(i) = 1;
end

%% Parameters
MAX_TOL = 1e-10;
ROUND_COUNT_ON = 0;

%% Main Program
count = 0;
current_CHAIN = [];
if ROUND_COUNT_ON == 1
    round_count = zeros(2^N, 1);
end

for i=1:N
    next_arrival_time(i) = get_interarrival_time(arrival_rate(i),number_of_mini_slots_per_unit_time);
end
while currentT < simT
    count = count + 1;
    round_idx = 1;
    if count > 1
        total_delivery_history(:,count) = total_delivery_history(:,count-1);
    end
    if count >= Npoints
        qn_history = [qn_history zeros(N, Nmorepoints)];
        total_delivery_history = [total_delivery_history zeros(N, Nmorepoints)];
        frame_timestamps = [frame_timestamps zeros(1, Nmorepoints)];
        Npoints = Npoints + Nmorepoints;
    end
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
        
        %% Determine contention time (including possible collisions)
        contention_time = DIFS;   
        % Each node is active if it has at least one packet is active if
        contention_set = find(qn >= 1);      
        % There can be multiple winners, a.k.a. collision
        if ~isempty(contention_set)
            total_contention_interval = total_contention_interval + 1;
            winner_id = contention_set(find(min(backoff(contention_set)) == backoff(contention_set)));
            contention_time = contention_time + min(backoff(contention_set))*slot_time;
            backoff(contention_set) = backoff(contention_set) - min(backoff(contention_set));
            if length(winner_id) < 1 + MAX_TOL
                WIFI_nodes{winner_id(1)}.CW = WIFI_nodes{winner_id(1)}.CW_Min;
                backoff(winner_id(1)) = get_backoff(WIFI_nodes{winner_id(1)}.CW); 
            else
                for j=1:length(winner_id)
                    WIFI_nodes{winner_id(j)}.CW = min((WIFI_nodes{winner_id(j)}.CW)*2, WIFI_nodes{winner_id(j)}.CW_Max);
                    backoff(winner_id(j)) = get_backoff(WIFI_nodes{winner_id(j)}.CW); 
                end
                contention_time = contention_time + Tpkt; % time wasted during collision
                total_collision = total_collision + 1;
                winner_id = -2; % no winner due to collision
            end
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
                round_idx = round_idx + power(2,winner_id-1);
                % channel is unreliable
                if rand(1) <= channel_pn(winner_id)
                    qn(winner_id) = max(0,qn(winner_id) - 1); 
                    total_delivery_history(winner_id,count) = total_delivery_history(winner_id,count) + 1;
                    Tput_shortterm(winner_id) = Tput_shortterm(winner_id) + 1;                 
                    next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                    while next_in_CHAIN ~= winner_id
                        if qn(next_in_CHAIN) > 0
                            currentT = currentT + Tpkt; 
                            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                            time_since_last_queue_update = currentT;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                            % channel is unreliable
                            if rand(1) <= channel_pn(next_in_CHAIN)
                                qn(next_in_CHAIN) = max(0,qn(next_in_CHAIN) - 1);
                                total_delivery_history(next_in_CHAIN,count) = total_delivery_history(next_in_CHAIN,count) + 1;
                                Tput_shortterm(next_in_CHAIN) = Tput_shortterm(next_in_CHAIN) + 1;
                                next_in_CHAIN = WIFI_nodes{next_in_CHAIN}.succ_id;
                            else
                                break
                            end
                        else
                            % CHAIN breaks and the next contention interval
                            % starts
                            break
                        end
                    end
                end
            else
                % The winner follows DCF
                currentT = currentT + Tpkt; 
                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                time_since_last_queue_update = currentT;
                round_idx = round_idx + power(2,winner_id-1); 
                % channel is unreliable
                if rand(1) <= channel_pn(next_in_CHAIN)                
                    qn(winner_id) = max(0,qn(winner_id) - 1);  % service is deterministic
                    total_delivery_history(winner_id,count) = total_delivery_history(winner_id,count) + 1;
                    Tput_shortterm(winner_id) = Tput_shortterm(winner_id) + 1; 
                end
            end
        else
            if winner_id == -2
                % Do nothing
            else
                if winner_id == -1
                    % there is no packet available in the network
                    currentT = min(next_arrival_time);
                    queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                    time_since_last_queue_update = currentT;
                end
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

        %% Determine contention time (including possible collisions)
        contention_time = DIFS;   
        % Each node is active if it has at least one packet is active if
        contention_set = find(qn >= 1);      
        % There can be multiple winners, a.k.a. collision
        if ~isempty(contention_set)
            total_contention_interval = total_contention_interval + 1;
            winner_id = contention_set(find(min(backoff(contention_set)) == backoff(contention_set)));
            contention_time = contention_time + min(backoff(contention_set))*slot_time;
            backoff(contention_set) = backoff(contention_set) - min(backoff(contention_set));
            if length(winner_id) < 1 + MAX_TOL
                WIFI_nodes{winner_id(1)}.CW = WIFI_nodes{winner_id(1)}.CW_Min;
                backoff(winner_id(1)) = get_backoff(WIFI_nodes{winner_id(1)}.CW); 
            else
                for j=1:length(winner_id)
                    WIFI_nodes{winner_id(j)}.CW = min((WIFI_nodes{winner_id(j)}.CW)*2, WIFI_nodes{winner_id(j)}.CW_Max);
                    backoff(winner_id(j)) = get_backoff(WIFI_nodes{winner_id(j)}.CW); 
                end
                contention_time = contention_time + Tpkt; % time wasted during collision
                total_collision = total_collision + 1;
                winner_id = -2; % no winner due to collision
            end
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
                round_idx = round_idx + power(2,winner_id-1); 
                % channel is unreliable
                if rand(1) <= channel_pn(winner_id) 
                    qn(winner_id) = max(0,qn(winner_id) - 1); 
                    total_delivery_history(winner_id,count) = total_delivery_history(winner_id,count) + 1;
                    Tput_shortterm(winner_id) = Tput_shortterm(winner_id) + 1;
                    next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                    while next_in_CHAIN ~= winner_id
                        if qn(next_in_CHAIN) > 0
                            currentT = currentT + Tpkt; 
                            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                            time_since_last_queue_update = currentT;
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                            % channel is unreliable
                            if rand(1) <= channel_pn(next_in_CHAIN)
                                qn(next_in_CHAIN) = max(0,qn(next_in_CHAIN) - 1);
                                total_delivery_history(next_in_CHAIN,count) = total_delivery_history(next_in_CHAIN,count) + 1;
                                Tput_shortterm(next_in_CHAIN) = Tput_shortterm(next_in_CHAIN) + 1;
                                next_in_CHAIN = WIFI_nodes{next_in_CHAIN}.succ_id;
                            else
                                break
                            end
                        else
                            % CHAIN breaks and the next contention interval
                            % starts
                            break
                        end
                    end
                end
            else
                % The winner follows DCF
                currentT = currentT + Tpkt; 
                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                time_since_last_queue_update = currentT;
                round_idx = round_idx + power(2,winner_id-1); 
                % channel is unreliable
                if rand(1) <= channel_pn(winner_id) 
                    qn(winner_id) = max(0,qn(winner_id) - 1);
                    total_delivery_history(winner_id,count) = total_delivery_history(winner_id,count) + 1;
                    Tput_shortterm(winner_id) = Tput_shortterm(winner_id) + 1;   
                end
            end
        else
            if winner_id == -2
                % Do nothing
            else
                if winner_id == -1
                    % there is no packet available in the network
                    currentT = min(next_arrival_time);
                    queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                    time_since_last_queue_update = currentT;
                end
            end
        end
    end
    
%% Q-threshold plus contention-based re-admission
% Chain may break due to a link with an empty queue
    if strcmp(Mode,'Qth-plus-Contention')
        % Contention
        winner_id = -1;
        CHAIN_break = 0;
        %% Determine contention time (including possible collisions)
        contention_time = DIFS;   
        % CHAIN node is active if
        % (1) it is already in CHAIN and has at least one packet
        % (2) it it currently not in CHAIN but has more than or equal to
        % q_threshold packets
        contention_set_part1 = find(qn(1:N_CHAIN_node).*(state_vec(1:N_CHAIN_node) == 0) >= q_threshold);
        contention_set_part2 = find(qn(1:N_CHAIN_node).*(state_vec(1:N_CHAIN_node) == 1) >= 1);
        % DCF node is active if it has at least one packet or has saturated
        % traffic
        contention_set_part3 = find((qn.*is_DCF) + is_DCFsat >= 1);
        contention_set = [contention_set_part1; contention_set_part2; contention_set_part3];     
        % There can be multiple winners, a.k.a. collision
        if ~isempty(contention_set)
            total_contention_interval = total_contention_interval + 1;
            winner_id = contention_set(find(min(backoff(contention_set)) == backoff(contention_set)));
            contention_time = contention_time + min(backoff(contention_set))*slot_time;
            backoff(contention_set) = backoff(contention_set) - min(backoff(contention_set));
            if length(winner_id) < 1 + MAX_TOL
                WIFI_nodes{winner_id(1)}.CW = WIFI_nodes{winner_id(1)}.CW_Min;
                backoff(winner_id(1)) = get_backoff(WIFI_nodes{winner_id(1)}.CW); 
            else
                for j=1:length(winner_id)
                    WIFI_nodes{winner_id(j)}.CW = min((WIFI_nodes{winner_id(j)}.CW)*2, WIFI_nodes{winner_id(j)}.CW_Max);
                    backoff(winner_id(j)) = get_backoff(WIFI_nodes{winner_id(j)}.CW); 
                end
                contention_time = contention_time + Tpkt; % time wasted during collision
                total_collision = total_collision + 1;
                winner_id = -2; % no winner due to collision
            end
            currentT = currentT + contention_time;
            total_contention_time = total_contention_time + contention_time;
            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
            time_since_last_queue_update = currentT;
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
                    round_idx = round_idx + power(2,winner_id-1);
                    % channel is unreliable
                    if rand(1) <= channel_pn(winner_id) 
                        qn(winner_id) = max(0,qn(winner_id) - 1);
                        total_delivery_history(winner_id,count) = total_delivery_history(winner_id,count) + 1;
                        next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                        while next_in_CHAIN ~= winner_id
                            if qn(next_in_CHAIN) > 0
                                currentT = currentT + Tpkt; 
                                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                                time_since_last_queue_update = currentT;
                                round_idx = round_idx + power(2,next_in_CHAIN-1);
                                % channel is unreliable
                                if rand(1) <= channel_pn(next_in_CHAIN) 
                                    qn(next_in_CHAIN) = max(0,qn(next_in_CHAIN) - 1);
                                    total_delivery_history(next_in_CHAIN,count) = total_delivery_history(next_in_CHAIN,count) + 1;
                                else
                                    CHAIN_break = 1;
                                end
                            else
                                CHAIN_break = 1;
                            end
                            % CHAIN breaks and needs to update predecessors
                            % and successors
                            if CHAIN_break == 1
                                WIFI_nodes{next_in_CHAIN}.is_active = 0;
                                state_vec(next_in_CHAIN) = 0;
                                current_CHAIN = current_CHAIN(current_CHAIN ~= next_in_CHAIN);
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
                        % No need to change Q-CHAIN relations                       
                    end
                else
                    % The winner is currently not in CHAIN
                    % TODO: The winner joins the CHAIN as the head, and the
                    % rest of CHAIN can piggyback
                    currentT = currentT + Tpkt; 
                    queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                    time_since_last_queue_update = currentT;
                    round_idx = round_idx + power(2,winner_id-1);
                    % channel is unreliable
                    if rand(1) <= channel_pn(winner_id) 
                        WIFI_nodes{winner_id}.is_active = 1;
                        current_CHAIN = [winner_id current_CHAIN];
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
                        qn(winner_id) = max(0,qn(winner_id) - 1); 
                        total_delivery_history(winner_id,count) = total_delivery_history(winner_id,count) + 1;
                    
                        % The rest of CHAIN can piggyback with possibly dummy
                        % packets
                        next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                        while next_in_CHAIN ~= winner_id
                            if qn(next_in_CHAIN) > 0
                                currentT = currentT + Tpkt; 
                                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                                time_since_last_queue_update = currentT;
                                round_idx = round_idx + power(2,next_in_CHAIN-1);
                                % channel is unreliable
                                if rand(1) <= channel_pn(next_in_CHAIN)
                                    qn(next_in_CHAIN) = max(0,qn(next_in_CHAIN) - 1);
                                    total_delivery_history(next_in_CHAIN,count) = total_delivery_history(next_in_CHAIN,count) + 1;
                                else
                                    % When a new node joins and the CHAIN
                                    % breaks in the middle, all the nodes in the
                                    % remaining part of the CHAIN should
                                    % leave
                                    % Remove one by one
                                    while next_in_CHAIN ~= CHAIN_head_id
                                        WIFI_nodes{next_in_CHAIN}.is_active = 0;
                                        state_vec(next_in_CHAIN) = 0;
                                        current_CHAIN = current_CHAIN(current_CHAIN ~= next_in_CHAIN);
                                        % Update CHAIN head or tail
                                        if next_in_CHAIN == CHAIN_tail_id
                                            CHAIN_tail_id = WIFI_nodes{next_in_CHAIN}.pred_id;
                                        end
                                        next_in_CHAIN_new = WIFI_nodes{next_in_CHAIN}.succ_id;
                                        WIFI_nodes{WIFI_nodes{next_in_CHAIN}.succ_id}.pred_id = WIFI_nodes{next_in_CHAIN}.pred_id;
                                        WIFI_nodes{WIFI_nodes{next_in_CHAIN}.pred_id}.succ_id = WIFI_nodes{next_in_CHAIN}.succ_id;
                                        WIFI_nodes{next_in_CHAIN}.succ_id = 0;
                                        WIFI_nodes{next_in_CHAIN}.pred_id = 0;
                                        next_in_CHAIN = next_in_CHAIN_new;
                                    end
                                    break
                                end
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
                    else
                        % Do nothing
                    end
                end
            else
                % The winner follows DCF
                currentT = currentT + Tpkt; 
                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                time_since_last_queue_update = currentT;
                round_idx = round_idx + power(2,winner_id-1);
                % channel is unreliable
                if rand(1) <= channel_pn(winner_id) 
                    qn(winner_id) = max(0,qn(winner_id) - 1);
                    total_delivery_history(winner_id,count) = total_delivery_history(winner_id,count) + 1;
                end
            end 
        else
            if winner_id == -2
                % Do nothing
            else
                if winner_id == -1
                    % there is no packet available in the network
                    currentT = min(next_arrival_time);
                    queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                    time_since_last_queue_update = currentT;
                end
            end
        end
    end

%% Q-threshold plus contention-based re-admission
% Chain may break due to a link with an empty queue
% A link can stay in the Chain by sending dummy packets
% NOTE: Unreliable channel for this policy has not been implemented yet!!!!
    if strcmp(Mode,'Qth-plus-Contention-with-dummy')
        % Contention
        winner_id = -1;

        %% Determine contention time (including possible collisions)
        contention_time = DIFS;   
        % CHAIN node is active if
        % (1) it is already in CHAIN and has at least one packet
        % (2) it it currently not in CHAIN but has more than or equal to
        % q_threshold packets
        contention_set_part1 = find(qn(1:N_CHAIN_node).*(state_vec(1:N_CHAIN_node) == 0) >= q_threshold);
        contention_set_part2 = find(qn(1:N_CHAIN_node).*(state_vec(1:N_CHAIN_node) == 1) >= 1);
        % DCF node is active if it has at least one packet
        contention_set_part3 = find(qn.*is_DCF >= 1);
        contention_set = [contention_set_part1; contention_set_part2; contention_set_part3];     
        % There can be multiple winners, a.k.a. collision
        if ~isempty(contention_set)
            total_contention_interval = total_contention_interval + 1;
            winner_id = contention_set(find(min(backoff(contention_set)) == backoff(contention_set)));
            contention_time = contention_time + min(backoff(contention_set))*slot_time;
            backoff(contention_set) = backoff(contention_set) - min(backoff(contention_set));
            if length(winner_id) < 1 + MAX_TOL
                WIFI_nodes{winner_id(1)}.CW = WIFI_nodes{winner_id(1)}.CW_Min;
                backoff(winner_id(1)) = get_backoff(WIFI_nodes{winner_id(1)}.CW); 
            else
                for j=1:length(winner_id)
                    WIFI_nodes{winner_id(j)}.CW = min((WIFI_nodes{winner_id(j)}.CW)*2, WIFI_nodes{winner_id(j)}.CW_Max);
                    backoff(winner_id(j)) = get_backoff(WIFI_nodes{winner_id(j)}.CW); 
                end
                contention_time = contention_time + Tpkt; % time wasted during collision
                total_collision = total_collision + 1;
                winner_id = -2; % no winner due to collision
            end
            currentT = currentT + contention_time;
            total_contention_time = total_contention_time + contention_time;
            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
            time_since_last_queue_update = currentT;
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
                    qn(winner_id) = max(0,qn(winner_id) - 1);  % service is deterministic
                    round_idx = round_idx + power(2,winner_id-1);
                    next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                    while next_in_CHAIN ~= winner_id
                        if qn(next_in_CHAIN) > 0
                            currentT = currentT + Tpkt; 
                            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                            time_since_last_queue_update = currentT;
                            qn(next_in_CHAIN) = max(0,qn(next_in_CHAIN) - 1);  % service is deterministic
                            round_idx = round_idx + power(2,next_in_CHAIN-1);
                        else
                            % CHAIN breaks and needs to update predecessors
                            % and successors
                            if dummy_packet_count(next_in_CHAIN) >= dummy_packet_limit
                                dummy_packet_count(next_in_CHAIN) = 0; % reset dummy packet counter
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
                            else
                                % Send a dummy packet to stay in Chain
                                currentT = currentT + Tdummy;
                                dummy_packet_count(next_in_CHAIN) = dummy_packet_count(next_in_CHAIN) + 1;
                                queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                                time_since_last_queue_update = currentT;
                                round_idx = round_idx + power(2,next_in_CHAIN-1);
                            end
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
                    qn(winner_id) = max(0,qn(winner_id) - 1);  % service is deterministic
                    round_idx = round_idx + power(2,winner_id-1);
                    
                    % The rest of CHAIN can piggyback with possibly dummy
                    % packets
                    next_in_CHAIN = WIFI_nodes{winner_id}.succ_id;
                    while next_in_CHAIN ~= winner_id
                        if qn(next_in_CHAIN) > 0
                            currentT = currentT + Tpkt; 
                            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                            time_since_last_queue_update = currentT;
                            qn(next_in_CHAIN) = max(0,qn(next_in_CHAIN) - 1);  % service is deterministic
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
                qn(winner_id) = max(0,qn(winner_id) - 1);  % service is deterministic
                round_idx = round_idx + power(2,winner_id-1);
            end 
        else
            if winner_id == -2
                % Do nothing
            else
                if winner_id == -1
                    % there is no packet available in the network
                    currentT = min(next_arrival_time);
                    queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                    time_since_last_queue_update = currentT;
                end
            end
        end
    end    

%% Queue-Length-Based CSMA (TON 2011, Jiang and Walrand)
    if strcmp(Mode,'QLB-CSMA')
        % Contention
        winner_id = -1;
        
        % Update queue length periodically
        if (currentT - time_last_queue_update_QLBCSMA) >= qlen_update_interval
            qn_slow = qn;
            time_last_queue_update_QLBCSMA = currentT;
        end
        %% Determine contention time (including possible collisions)
        contention_time = DIFS;   
        % Each node is active if it has at least one packet is active if
        contention_set = find((qn >= 1) | (is_DCFsat == 1));      
        % There can be multiple winners, a.k.a. collision
        if ~isempty(contention_set)
            total_contention_interval = total_contention_interval + 1;
            winner_id = contention_set(find(min(backoff(contention_set)) == backoff(contention_set)));
            contention_time = contention_time + min(backoff(contention_set))*slot_time;
            backoff(contention_set) = backoff(contention_set) - min(backoff(contention_set));
            if length(winner_id) < 1 + MAX_TOL
                WIFI_nodes{winner_id(1)}.CW = WIFI_nodes{winner_id(1)}.CW_Min;
                backoff(winner_id(1)) = get_backoff(getCW_QLB(qn_slow(winner_id(1)))); 
            else
                for j=1:length(winner_id)
                    WIFI_nodes{winner_id(j)}.CW = min((WIFI_nodes{winner_id(j)}.CW)*2, WIFI_nodes{winner_id(j)}.CW_Max);
                    backoff(winner_id(j)) = get_backoff(getCW_QLB(qn_slow(winner_id(j)))); 
                end
                contention_time = contention_time + Tpkt; % time wasted during collision
                total_collision = total_collision + 1;
                winner_id = -2; % no winner due to collision
            end
            currentT = currentT + contention_time;
            total_contention_time = total_contention_time + contention_time;
            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
            time_since_last_queue_update = currentT;
        end
        if winner_id > 0            
            % Time evolution
            currentT = currentT + Tpkt; 
            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
            time_since_last_queue_update = currentT;
            round_idx = round_idx + power(2,winner_id-1);
            % channel is unreliable
            if rand(1) <= channel_pn(winner_id)
                qn(winner_id) = max(0,qn(winner_id) - 1); 
                total_delivery_history(winner_id,count) = total_delivery_history(winner_id,count) + 1;
            end
        else
            if winner_id == -2
                % Do nothing
            else
                if winner_id == -1
                    % there is no packet available in the network
                    currentT = min(next_arrival_time);
                    queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                    time_since_last_queue_update = currentT;
                end
            end           
        end
    end
%% Conventional Pure-DCF mode
    if strcmp(Mode,'Pure-DCF')
        % Contention
        winner_id = -1;
        
        %% Determine contention time (including possible collisions)
        contention_time = DIFS;   
        % Each node is active if it has at least one packet is active if
        contention_set = find((qn >= 1) | (is_DCFsat == 1));      
        % There can be multiple winners, a.k.a. collision
        if ~isempty(contention_set)
            total_contention_interval = total_contention_interval + 1;
            winner_id = contention_set(find(min(backoff(contention_set)) == backoff(contention_set)));
            contention_time = contention_time + min(backoff(contention_set))*slot_time;
            backoff(contention_set) = backoff(contention_set) - min(backoff(contention_set));
            if length(winner_id) < 1 + MAX_TOL
                WIFI_nodes{winner_id(1)}.CW = WIFI_nodes{winner_id(1)}.CW_Min;
                backoff(winner_id(1)) = get_backoff(WIFI_nodes{winner_id(1)}.CW); 
            else
                for j=1:length(winner_id)
                    WIFI_nodes{winner_id(j)}.CW = min((WIFI_nodes{winner_id(j)}.CW)*2, WIFI_nodes{winner_id(j)}.CW_Max);
                    backoff(winner_id(j)) = get_backoff(WIFI_nodes{winner_id(j)}.CW); 
                end
                contention_time = contention_time + Tpkt; % time wasted during collision
                total_collision = total_collision + 1;
                winner_id = -2; % no winner due to collision
            end
            currentT = currentT + contention_time;
            total_contention_time = total_contention_time + contention_time;
            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
            time_since_last_queue_update = currentT;
        end            
        
        if winner_id > 0            
            % Time evolution
            currentT = currentT + Tpkt; 
            queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
            time_since_last_queue_update = currentT;
            round_idx = round_idx + power(2,winner_id-1);
            % channel is unreliable
            if rand(1) <= channel_pn(winner_id)
                qn(winner_id) = max(0,qn(winner_id) - 1); 
                total_delivery_history(winner_id,count) = total_delivery_history(winner_id,count) + 1;
            end
        else
            if winner_id == -2
                % Do nothing
            else
                if winner_id == -1
                    % there is no packet available in the network
                    currentT = min(next_arrival_time);
                    queue_length_calculus = queue_length_calculus + (currentT - time_since_last_queue_update)*qn;
                    time_since_last_queue_update = currentT;
                end
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
    if ROUND_COUNT_ON == 1
        round_count(round_idx,1) = round_count(round_idx,1) + 1;
    end
    % Update history
    qn_history(:,count) = qn;
    frame_timestamps(count) = currentT;
end

%% Post-processing
if ROUND_COUNT_ON == 1
    round_count_normalize = round_count/(sum(round_count));
end
average_contention_time_per_interval = total_contention_time/total_contention_interval;
average_throughput = max(total_delivery_history, [], 2)/currentT;

%% Plotting
%plot(frame_timestamps(1:count), qn_history(:,1:count));
Stepsize = 50;
createfigure(frame_timestamps(1:Stepsize:count), qn_history(:,1:Stepsize:count));
average_qn = queue_length_calculus/currentT;
switch N_arrival_groups
    case 1
        fprintf(fileID, 'average throughput of Q-CHAIN clients (pkt/sec): %.2f\n', sum(average_throughput(1:N_CHAIN_node))*1000);
        fprintf(fileID, 'average throughput of DCF clients (pkt/sec): %.2f\n', sum(average_throughput(N_CHAIN_node+1:N_CHAIN_node+N_DCF))*1000);
        fprintf(fileID, 'average throughput of DCFsat clients (pkt/sec): %.2f\n', sum(average_throughput(N_CHAIN_node+N_DCF+1:N_CHAIN_node+N_DCF+N_DCFsat))*1000);
        fprintf(fileID, 'total average queue length: %.2f\n', sum(average_qn));
        fprintf(fileID, 'total average delay: %.2f (ms)\n', sum(average_qn)/sum(arrival_rate));
    case 2
        fprintf(fileID, 'average throughput of Q-CHAIN clients (pkt/sec): %.2f\n', sum(average_throughput(1:N_CHAIN_node))*1000);
        fprintf(fileID, 'average throughput of DCF clients (pkt/sec): %.2f\n', sum(average_throughput(N_CHAIN_node+1:N_CHAIN_node+N_DCF))*1000);
        fprintf(fileID, 'average throughput of DCFsat clients (pkt/sec): %.2f\n', sum(average_throughput(N_CHAIN_node+N_DCF+1:N_CHAIN_node+N_DCF+N_DCFsat))*1000);
        fprintf(fileID, 'average throughput of Arrival-Group-1 clients (pkt/sec): %.2f\n', sum(average_throughput(1:N_half))*1000);
        fprintf(fileID, 'average throughput of Arrival-Group-2 clients (pkt/sec): %.2f\n', sum(average_throughput(N_half+1:N_nonsat))*1000);
        fprintf(fileID, 'total average queue length of Arrival-Group-1: %.2f\n', sum(average_qn(1:N_half)));
        fprintf(fileID, 'total average queue length of Arrival-Group-2: %.2f\n', sum(average_qn(N_half+1:N_nonsat)));
        fprintf(fileID, 'total average queue length: %.2f\n', sum(average_qn));
        fprintf(fileID, 'total average delay of Arrival-Group-1: %.2f (ms)\n', sum(average_qn(1:N_half))/sum(arrival_rate(1:N_half)));
        fprintf(fileID, 'total average delay of Arrival-Group-2: %.2f (ms)\n', sum(average_qn(N_half+1:N_nonsat))/sum(arrival_rate(N_half+1:N_nonsat)));
        fprintf(fileID, 'total average delay: %.2f (ms)\n', sum(average_qn)/sum(arrival_rate));
    otherwise
        fprintf(fileID, 'average throughput of Q-CHAIN clients (pkt/sec): %.2f\n', sum(average_throughput(1:N_CHAIN_node))*1000);
        fprintf(fileID, 'average throughput of DCF clients (pkt/sec): %.2f\n', sum(average_throughput(N_CHAIN_node+1:N_CHAIN_node+N_DCF))*1000);
        fprintf(fileID, 'average throughput of DCFsat clients (pkt/sec): %.2f\n', sum(average_throughput(N_CHAIN_node+N_DCF+1:N_CHAIN_node+N_DCF+N_DCFsat))*1000);
        fprintf(fileID, 'total average queue length: %.2f\n', sum(average_qn));
        fprintf(fileID, 'total average delay: %.2f (ms)\n', sum(average_qn)/sum(arrival_rate));
end
fprintf(fileID, 'average contention time per interval: %.3f (ms)\n', average_contention_time_per_interval);
fprintf(fileID, 'contention time percentage: %.3f %%\n', total_contention_time/currentT*100);
fprintf(fileID, 'average collision rate: %.2f %%\n\n', total_collision/total_contention_interval*100);
for i=1:N
fprintf(fileID, 'average queue length of client %d: %.2f\n', i, average_qn(i));
end
for i=1:N
fprintf(fileID, 'average throughput of client %d (pkt/sec): %.2f\n', i, average_throughput(i)*1000);
end


toc;
end

