%% Plots for average delay and contention time
clear;
clc;
filepath_delay = '../figures/delay/';
filepath_contention_percentage = '../figures/contention-percentage/';
filepath_contention_time_per_interval = '../figures/contention-time-per-interval/';
filepath_collision_percentage = '../figures/collision-percentage/';

case_code = '20-1-10-2groups-het';

switch case_code
    
    case '20-1-0-2groups-het'
%% Case: 20 nodes, 1 Chain, 0 DCF nodes, 2 groups, heterogeneous arrival rates
%% Qth-plus-Contention
Qth_plus_Contention_rho = [0.7 0.72 0.74 0.76 0.78 0.8 0.81 0.82 0.83 0.84 0.85 0.86];
Qth_plus_Contention_total_delay_ms = [21.19 20.92 20.98 21.00 21.43 22.55 24.20 26.86 29.00 38.33 57.46 277.22];
Qth_plus_Contention_contention_per_interval_us = 1000*[0.089 0.089 0.088 0.088 0.088 0.087 0.088 0.088 0.088 0.087 0.088 0.088];
Qth_plus_Contention_contention_time_percentage = [19.90 19.24 18.39 17.67 16.56 15.37 14.58 13.78 12.96 11.8 11.05 10.08];
Qth_plus_Contention_collision_percentage = [8.32 9.37 10.40 11.45 12.97 14.31 15.56 16.64 17.56 18.79 19.97 20.99];



%% Original CHAIN-1
CHAIN_1_rho = [0.6 0.62 0.64 0.65 0.66 0.67 0.68 0.69 0.70 0.71];
CHAIN_1_total_delay_ms = [1.36 1.81 2.54 3.25 4.24 7.50 12.30 22.58 35.50 426.14];
CHAIN_1_contention_per_interval_us = 1000*[0.089 0.089 0.089 0.088 0.088 0.087 0.087 0.087 0.087 0.087];
CHAIN_1_contention_time_percentage = [30.97 31.02 30.73 30.36 30.03 29.39 28.78 27.87 27.08 26.32];
CHAIN_1_collision_percentage = [10.03 11.60 13.22 14.10 14.99 16.09 16.76 17.66 18.26 18.93];



%% Original CHAIN-2
CHAIN_2_rho = [0.6 0.62 0.64 0.65 0.66 0.67 0.68 0.69 0.70 0.71];
CHAIN_2_total_delay_ms = [1.37 1.76 2.61 3.12 4.81 7.34 14.95 23.92 45.39 349.49];
CHAIN_2_contention_per_interval_us = 1000*[0.089 0.088 0.089 0.088 0.088 0.088 0.087 0.087 0.087 0.087];
CHAIN_2_contention_time_percentage = [30.99 30.98 30.77 30.55 30.12 29.52 28.76 27.95 26.85 26.29];
CHAIN_2_collision_percentage = [10.14 11.47 13.44 14.23 15.00 15.97 16.97 17.54 18.45 18.98];



%% DCF
DCF_rho = [0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58];
DCF_total_delay_ms = [0.95 1.13 1.32 1.55 2.13 3.33 6.22 24.25 933.02];
DCF_contention_per_interval_us = 1000*[0.091 0.090 0.090 0.089 0.089 0.089 0.088 0.088 0.087];
DCF_contention_time_percentage = [31.86 32.71 33.43 34.19 35.31 36.55 37.96 39.32 40.43];
DCF_collision_percentage = [7.21 8.10 8.89 9.67 11.02 12.68 14.70 17.17 19.50];



    case '20-1-10-2groups-het'
%% Case: 20 nodes, 1 Chain, 10 DCF nodes, 2 groups, heterogeneous arrival rates
%% Qth-plus-Contention
Qth_plus_Contention_rho = [0.66 0.68 0.7 0.72 0.74 0.75 0.76 0.77 0.78];
Qth_plus_Contention_total_delay_ms = [12.00 12.15 12.87 14.39 20.35 29.95 39.89 83.37 465.45];
Qth_plus_Contention_contention_per_interval_us = 1000*[0.089 0.089 0.090 0.090 0.090 0.091 0.091 0.092 0.092];
Qth_plus_Contention_contention_time_percentage = [24.85 24.58 24.21 23.50 22.20 21.28 20.53 19.59 19.08];
Qth_plus_Contention_collision_percentage = [10.90 12.31 14.26 16.20 19.14 21.31 22.37 24.31 25.03];


%% Original CHAIN-1
CHAIN_1_rho = [0.58 0.6 0.62 0.64 0.65 0.66 0.67 0.68 0.69];
CHAIN_1_total_delay_ms = [1.18 1.48 2.08 3.12 4.19 6.28 15.19 23.57 154.04];
CHAIN_1_contention_per_interval_us = 1000*[0.090 0.089 0.089 0.089 0.088 0.088 0.088 0.088 0.088];
CHAIN_1_contention_time_percentage = [31.16 31.42 31.49 31.17 30.94 30.46 29.75 28.97 28.11];
CHAIN_1_collision_percentage = [9.29 10.58 12.35 14.23 14.89 16.11 17.30 18.03 19.02];


%% Original CHAIN-2
CHAIN_2_rho = [0.58 0.6 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69];
CHAIN_2_total_delay_ms = [1.19 1.52 2.12 2.52 3.05 4.37 8.22 11.28 22.98 158.96];
CHAIN_2_contention_per_interval_us = 1000*[0.089 0.089 0.089 0.089 0.088 0.088 0.088 0.088 0.088 0.088];
CHAIN_2_contention_time_percentage = [31.13 31.34 31.52 31.34 31.12 30.89 30.34 29.73 28.99 28.10];
CHAIN_2_collision_percentage = [9.25 10.71 12.42 13.15 14.11 15.25 16.37 17.26 18.06 19.39];


%% DCF
DCF_rho = [0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58];
DCF_total_delay_ms = [0.95 1.13 1.32 1.55 2.13 3.33 6.22 24.25 933.02];
DCF_contention_per_interval_us = 1000*[0.091 0.090 0.090 0.089 0.089 0.089 0.088 0.088 0.087];
DCF_contention_time_percentage = [31.86 32.71 33.43 34.19 35.31 36.55 37.96 39.32 40.43];
DCF_collision_percentage = [7.21 8.10 8.89 9.67 11.02 12.68 14.70 17.17 19.50];

end


%%%%% Plotting
createfigure_delay(Qth_plus_Contention_rho, Qth_plus_Contention_total_delay_ms, ...
    CHAIN_1_rho, CHAIN_1_total_delay_ms, CHAIN_2_rho, CHAIN_2_total_delay_ms, ...
    DCF_rho, DCF_total_delay_ms, case_code, filepath_delay);
createfigure_contention_percentage(Qth_plus_Contention_rho, Qth_plus_Contention_contention_time_percentage, ...
    CHAIN_1_rho, CHAIN_1_contention_time_percentage, CHAIN_2_rho, CHAIN_2_contention_time_percentage, ...
    DCF_rho, DCF_contention_time_percentage, case_code, filepath_contention_percentage);
createfigure_contention_time_per_interval(Qth_plus_Contention_rho, Qth_plus_Contention_contention_per_interval_us, ...
    CHAIN_1_rho, CHAIN_1_contention_per_interval_us, CHAIN_2_rho, CHAIN_2_contention_per_interval_us, ...
    DCF_rho, DCF_contention_per_interval_us, case_code, filepath_contention_time_per_interval);
createfigure_collision_percentage(Qth_plus_Contention_rho, Qth_plus_Contention_collision_percentage, ...
    CHAIN_1_rho, CHAIN_1_collision_percentage, CHAIN_2_rho, CHAIN_2_collision_percentage, ...
    DCF_rho, DCF_collision_percentage, case_code, filepath_collision_percentage);


