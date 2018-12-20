%%% Main program for multiple simulation runs
clc; clear;
Mode = 'Qth-plus-Contention';
%{'Pure-DCF', 'Original-CHAIN', 'Original-CHAIN-two-chains', 'Qth-plus-Contention'};
%rho_vec = [0.7 0.72 0.74 0.76 0.78 0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.93 0.94];
%rho_vec = [0.008];
rho_vec = [0.013 0.012 0.011 0.01 0.009 0.008 0.007 0.006 0.005 0.004 0.003 0.002 0.001];
%rho_vec = [0.007 0.006 0.005 0.004 0.003 0.002 0.001];
%rho_vec = 0.013;
%rho_vec = 0.65;
%network_config_code = '2-1-0-0-1group';
%network_config_code = '3-1-0-0-1group';
%network_config_code = '10-1-0-0-1group';
%network_config_code = '20-1-10-2groups-het';
network_config_code = '500-1-0-0-1group';
%network_config_code = '500-1-0-0-2groups';

for i=1:length(rho_vec)
    pause(1)
    time = datestr(datetime('now','TimeZone','local','Format','d-MMM-y_HH_mm_ss'));
    fileID = fopen(strcat('./rawdata/sim_data_', network_config_code, '/', network_config_code, '_', Mode, '_', time),'w');
    sim_single_CHAINs_real_backoff(network_config_code, Mode, fileID, rho_vec(i));
end