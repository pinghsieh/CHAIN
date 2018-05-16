%%% Main program for multiple simulation runs

Mode = 'Qth-plus-Contention';
%{'Pure-DCF', 'Original-CHAIN', 'Original-CHAIN-two-chains', 'Qth-plus-Contention'};
rho1_vec = [0.3];
rho_vec = [0.6];
network_config_code = '20-1-0-2groups-cap';
time = datestr(datetime('now','TimeZone','local','Format','d-MMM-y_HH_mm_ss'));
fileID = fopen(strcat('./rawdata/sim_data_', network_config_code, '_', Mode, '_', time),'w');

for i=1:length(rho1_vec)
    rho2_vec = rho_vec - rho1_vec(i);
    for j=1:length(rho2_vec)
        sim_single_CHAINs_real_backoff_capacity_region(network_config_code, rho1_vec(i), rho2_vec(j), Mode, fileID);
    end
end