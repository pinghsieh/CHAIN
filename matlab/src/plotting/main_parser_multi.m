%% Parser for end-to-end delay log
% real_delay: unit is millisecond
% deadline: unit is millisecond
% relative_delay: unit is 0.1 microsecond
N = 6;
filename = cell(N,1);
real_delay = cell(N,1);
cdf_delay = cell(N,1);
pdf_delay = cell(N,1);
Xvalue = cell(N,1);
sample = cell(N,1);
Num_ACK_entry = zeros(N,1);
Num_ACK_timeout = zeros(N,1);
mean_delay = zeros(N,1);
delimiterIn = ' ';
headerlinesIn = 0;
MIN_TOL = 1e-5;
spacing = 50;
%filename{1} = './Test Results 20170426/latency_500us.txt';
%filename{2} = './Test Results 20170426/latency_1ms.txt';
%filename{3} = './Test Results 20170426/latency_2ms.txt';
%filename{4} = './Test Results 20170426/latency_5ms.txt';
%filename{1} = './Test Results 20170426/lifetime_500us.txt';
%filename{2} = './Test Results 20170426/lifetime_1ms.txt';
%filename{3} = './Test Results 20170426/lifetime_2ms.txt';
%filename{4} = './Test Results 20170426/lifetime_5ms.txt';
%filename{1} = './Test Results 20170725/lifetime.txt';
%directory{1} = './Test Results 20170726/latency_50k_MCS7_1500B_10ms/';
%directory{1} = './Test Results 20170727 MultiStation/';
directory{1} = './Test Results 20170729/10ms_Periodic_Final/';
MCS_vec = [7 7 7 4 4 4];
Rate_vec = [54, 54, 54, 24, 24, 24];
Payload_vec = [500 1000 1500 500 1000 1500];
%mode = 'Host-to-FPGA';
mode = 'Round-Trip';
%trial = '5';
%filename{1} = strcat(directory{1}, 'latency_trial1.txt');
%filename{1} = strcat(directory{1}, 'lifetime_trial5.txt');
%filename{1} = strcat(directory{1}, 'latency_MCS4_1500B.txt');
%{
filename{1} = strcat(directory{1}, 'lifetime_MCS7_500B.txt');
filename{2} = strcat(directory{1}, 'lifetime_MCS7_1000B.txt');
filename{3} = strcat(directory{1}, 'lifetime_MCS7_1500B.txt');
filename{4} = strcat(directory{1}, 'lifetime_MCS4_500B.txt');
filename{5} = strcat(directory{1}, 'lifetime_MCS4_1000B.txt');
filename{6} = strcat(directory{1}, 'lifetime_MCS4_1500B.txt');
%}

filename{1} = strcat(directory{1}, 'latency_MCS7_500B.txt');
filename{2} = strcat(directory{1}, 'latency_MCS7_1000B.txt');
filename{3} = strcat(directory{1}, 'latency_MCS7_1500B.txt');
filename{4} = strcat(directory{1}, 'latency_MCS4_500B.txt');
filename{5} = strcat(directory{1}, 'latency_MCS4_1000B.txt');
filename{6} = strcat(directory{1}, 'latency_MCS4_1500B.txt');

%deadline = [0.5, 1, 2, 5];
deadline = [10 10 10 10 10 10]; %unit: 1ms
%color = {[0 0 1], [1 0 1], [1 0 0], [0 1 0]};
%color = {[0 0 1], [1 0 1]};
color = {[1 0 0], [1 0 1], [0 0 1], [0 1 0], [0 1 1], [1 0.75 0]};
marker = {'^', 's', 'o', 'x', 'd', '+'};
prc_vec = [90, 95, 99];

for i=1:N
    relative_delay = importdata(filename{i}, delimiterIn);
    K = length(relative_delay);
    ACK_entry = find(relative_delay > -2147483648 + MIN_TOL);
    Num_ACK_entry(i) = length(ACK_entry);
    Num_ACK_timeout(i) = K - length(ACK_entry);
    real_delay{i} = ((-1)*relative_delay(ACK_entry) + deadline(i)*10000)/10000;
    mean_delay(i) = mean(real_delay{i});
    delay_prc_value = prctile(real_delay{i}, prc_vec);
end

Num_sent = Num_ACK_entry + Num_ACK_timeout;

for i=1:N
   [cdf_delay{i},Xvalue{i}] = ecdf(real_delay{i});
   pdf_delay{i} = diff([0; cdf_delay{i}]);
   sample{i} = 1:spacing:length(pdf_delay{i});
   %createfigure(x,pdf_delay, deadline(i), color, mode, MCS, Payload);
   %createfigure(x,cdf_delay, deadline(i), color, mode, MCS, Payload);
   %savefig(strcat(directory{1}, mode, trial, '.fig'));
   %saveas(gcf, strcat(directory{1}, mode, trial, '.bmp'));
   %savefig(strcat(directory{1}, mode, '_MCS', num2str(MCS), '_', num2str(Payload), 'B','.fig'));
   %saveas(gcf, strcat(directory{1}, mode, '_MCS', num2str(MCS), '_', num2str(Payload), 'B', '.bmp'));
end

createfigure_multi(Xvalue, cdf_delay, sample, deadline, color, mode, Rate_vec, Payload_vec, marker, N);
savefig(strcat(directory{1}, mode, '_all','.fig'));
saveas(gcf, strcat(directory{1}, mode, '_all','.bmp'));