%%% Main program for multiple simulation runs

Mode = 'Original-CHAIN-two-chains';
%{'Pure-DCF', 'Original-CHAIN', 'Original-CHAIN-two-chains', 'Qth-plus-Contention'};
%rho_vec = [0.7 0.72 0.74 0.76 0.78 0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.93 0.94];
rho_vec = [0.69];
network_config_code = '20-1-10-2groups-het';
time = datestr(datetime('now','TimeZone','local','Format','d-MMM-y_HH_mm_ss'));
fileID = fopen(strcat('./rawdata/sim_data_', network_config_code, '_', Mode, '_', time),'w');

for i=1:length(rho_vec)
    sim_single_CHAINs_real_backoff(network_config_code, rho_vec(i), Mode, fileID);
end