%% Plots for average delay and contention time
%% For MobiHoc 2019 paper
clear;
clc;
filepath_delay = '../figures/mobihoc/delay/';
filepath_tput = '../figures/mobihoc/throughput/';
filepath_contention_percentage = '../figures/mobihoc/contention-percentage/';
filepath_collision_percentage = '../figures/mobihoc/collision-percentage/';

case_code = '500-1-0-1-1group-p=1';
LARGE_DELAY = 20; % for unstable cases

switch case_code
    
case '500-1-0-1-1group-p=1'
%% Case: 500 nodes, 1 Chain, 0 DCF nodes, 1 DCFsat node, 1 group
%% Qth-plus-Contention
Qth_plus_Contention_eta = [5, 4.5, 4, 3.5, 3, 2.5, 2, 1.5, 1, 0.5];
Qth_plus_Contention_DCFsat_tput = [1783.23, 1968.13, 2160.28, 2355.8, 2553.8, 2752.08, 2955.05, 3166.75, 3384.12, 3604.31];
Qth_plus_Contention_contention_time_percentage = [30.74, 31.77, 32.79, 33.77, 34.72, 35.55, 36.33, 37.09, 37.73, 38.34];


%% Original CHAIN-1
CHAIN_1_eta = [5, 4.5, 4, 3.5, 3, 2.5, 2, 1.5, 1, 0.5];
CHAIN_1_DCFsat_tput = [1319.84, 1631.51, 1910.66, 2176.07, 2427.79, 2664.92, 2910.98, 3143.94, 3372.36, 3603.19];
CHAIN_1_contention_time_percentage = [38.84, 37.96, 37.46, 37.19, 37.18, 37.3, 37.42, 37.72, 38.02, 38.37];

%% QLB-CSMA
QLB_CSMA_eta = [5, 4.5, 4, 3.5, 3, 2.5, 2, 1.5, 1, 0.5];
QLB_CSMA_DCFsat_tput = [352.13, 395.59, 434.78, 475.99, 513.38, 553.53, 591.85, 630.57, 669.23, 708.45];
QLB_CSMA_contention_time_percentage = [54.31, 57.8, 61.03, 64.44, 67.77, 71.18, 74.53, 77.93, 81.29, 84.67];

%% DCF
DCF_eta = [5, 4.5, 4, 3.5, 3, 2.5, 2, 1.5, 1, 0.5];
DCF_DCFsat_tput = [1297.93, 1624.11, 1904.43, 2171.44, 2421.3, 2667.18, 2907.14, 3141.04, 3374.4, 3600.34];
DCF_contention_time_percentage = [39.06, 38.06, 37.51, 37.24, 37.23, 37.27, 37.46, 37.71, 38.02, 38.38];

%%%%% Plotting

createfigure_tput_DCFsat_mobihoc(Qth_plus_Contention_eta, Qth_plus_Contention_DCFsat_tput, ...
    CHAIN_1_eta, CHAIN_1_DCFsat_tput, ...
    QLB_CSMA_eta, QLB_CSMA_DCFsat_tput, ...
    DCF_eta, DCF_DCFsat_tput, case_code, filepath_tput);

createfigure_contention_percentage_DCFsat_mobihoc(Qth_plus_Contention_eta, Qth_plus_Contention_contention_time_percentage, ...
    CHAIN_1_eta, CHAIN_1_contention_time_percentage, ...
    QLB_CSMA_eta, QLB_CSMA_contention_time_percentage, ...
    DCF_eta, DCF_contention_time_percentage, case_code, filepath_contention_percentage);

case '500-1-0-1-2groups-het-p=1'
%% Case: 500 nodes, 1 Chain, 0 DCF nodes, 1 DCFsat node, 2 groups, heterogeneous arrival rates, Kgroup1 = 2*Kgroup2
%% Qth-plus-Contention
Qth_plus_Contention_eta = [3.6, 3.5, 3, 2.5, 2, 1.5, 1, 0.5];
Qth_plus_Contention_DCFsat_tput = [1636.49, 1697.06, 1974.01, 2259.01 2551.27, 2851.02, 3163.37, 3493.47];
Qth_plus_Contention_contention_time_percentage = [29.831, 30.176, 31.769, 33.313, 34.708, 35.956, 37.089, 38.049];

%% Original CHAIN-1
CHAIN_1_eta = [3.6, 3.5, 3, 2.5, 2, 1.5, 1, 0.5];
CHAIN_1_DCFsat_tput = [1051.50, 1143.92, 1630.10, 2042.93, 2424.39, 2787.75, 3148.68, 3488.24];
CHAIN_1_contention_time_percentage = [40.045, 39.560, 37.905, 37.309, 37.206, 37.401, 37.702, 38.191];

%% QLB-CSMA
QLB_CSMA_eta = [3.6, 3.5, 3, 2.5, 2, 1.5, 1, 0.5];
QLB_CSMA_DCFsat_tput = [319.07, 332.61, 394.78, 454.41, 512.80, 574.41, 630.27, 687.52];
QLB_CSMA_contention_time_percentage = [51.748, 52.773, 57.747, 62.762, 67.750, 72.785, 77.840, 82.979];

%% DCF
DCF_eta = [3.6, 3.5, 3, 2.5, 2, 1.5, 1, 0.5];
DCF_DCFsat_tput = [1004.44, 1123.49, 1627.12, 2042.39, 2426.98, 2787.95, 3143.17, 3490.06];
DCF_contention_time_percentage = [40.721, 39.954, 37.971, 37.333, 37.195, 37.335, 37.718, 38.171];

%%%%% Plotting

createfigure_tput_DCFsat_mobihoc(Qth_plus_Contention_eta, Qth_plus_Contention_DCFsat_tput, ...
    CHAIN_1_eta, CHAIN_1_DCFsat_tput, ...
    QLB_CSMA_eta, QLB_CSMA_DCFsat_tput, ...
    DCF_eta, DCF_DCFsat_tput, case_code, filepath_tput);

createfigure_contention_percentage_DCFsat_mobihoc(Qth_plus_Contention_eta, Qth_plus_Contention_contention_time_percentage, ...
    CHAIN_1_eta, CHAIN_1_contention_time_percentage, ...
    QLB_CSMA_eta, QLB_CSMA_contention_time_percentage, ...
    DCF_eta, DCF_contention_time_percentage, case_code, filepath_contention_percentage);


case '500-1-0-0-2groups-het-p=0.9'
%% Case: 500 nodes, 1 Chain, 0 DCF nodes, 2 groups, heterogeneous arrival rates
%% Qth-plus-Contention
Qth_plus_Contention_rho = [6, 5.5, 5, 4.6, 4.5, 4, 3.5, 3, 2.5, 2];
Qth_plus_Contention_total_tput = [3920.26, 4002.46 4065.65, 4356.56, 4257.26, 3774.53, 3311.04, 2836.31, 2364.17, 1885.51];
Qth_plus_Contention_total_delay_s = [28.307, 21.136, 12.912, 0.642, 0.633, 0.67, 0.742, 0.847, 1.000, 1.230];
Qth_plus_Contention_contention_time_percentage = [29.91, 28.41, 27.24, 20.70, 21.27, 22.67, 22.49, 21.49, 19.64, 16.88];
Qth_plus_Contention_collision_percentage = [49.75, 47.25, 45.18, 18.55, 14.53, 8.13, 5.38, 3.95, 2.81, 1.95];



%% Original CHAIN-1
CHAIN_1_rho = [6, 5.5, 5, 4.5, 4, 3.5, 3.3, 3.2, 3, 2.5, 2];
CHAIN_1_total_tput = [3226.71, 3098.08, 3062.48, 3098.28, 3189.73, 3221.88, 3132.37, 3042.21, 2849.88, 2374.42, 1899.38];
CHAIN_1_total_delay_s = [39.012, 36.438, 31.986, 24.7, 14.463, 3.031, 0.392, 0.078, 2.2e-3, 0.56e-3, 0.41e-3];
CHAIN_1_contention_time_percentage = [42.61, 44.92, 45.55, 44.90, 43.31, 42.67, 44.28, 45.11, 33.94, 25.88, 20.72];
CHAIN_1_collision_percentage = [62.58, 61.34, 60.01, 58.84, 57.42, 54.19, 50.03, 40.30, 12.68, 4.65, 2.39];

%% QLB-CSMA
QLB_CSMA_rho = [6, 5.5, 5, 4.5, 4, 3.5, 3.3, 3.2, 3, 2.5, 2];
QLB_CSMA_total_tput = [1085.87, 1085.19, 1087.16 1087.21 1086.05, 1127.48, 1567.46, 3040.33, 2851.97, 2373.56, 1902.52];
QLB_CSMA_total_delay_s = [72.788, 71.176, 69.608, 67.023, 64.342, 58.354, 34.653, 5.74e-3, 4.66e-3, 3.28e-3, 2.62e-3];
QLB_CSMA_contention_time_percentage = [80.7, 80.71, 80.66, 80.69, 80.69, 79.97, 72.17, 45.94, 49.29, 57.70, 65.45];
QLB_CSMA_collision_percentage = [75.94, 75.95, 75.89, 75.92, 75.92, 74.94, 63.14, 6.07, 4.7, 2.79, 1.83];

%% DCF
DCF_rho = [6, 5.5, 5, 4.5, 4, 3.5, 3, 2.9, 2.8, 2.6, 2.5, 2];
DCF_total_tput = [1487.65, 1490.92, 1492.78, 1513.01, 1557.01, 1742.33, 2363.13, 2750.36, 2663.17, 2479.55, 2370.48, 1899.28];
DCF_total_delay_s = [66.433, 64.052, 61.543, 19.217,17.399, 13.239, 4.978, 1.67e-3, 0.99e-3, 0.63e-3, 0.58e-3, 0.41e-3];
DCF_contention_time_percentage = [73.55, 73.50, 73.48, 73.10, 72.37, 69.01, 57.19, 32.01, 29.85, 27.23, 26.01, 20.79];
DCF_collision_percentage = [66.27, 66.21, 66.18, 65.69, 64.75, 60.41, 45.19, 10.38, 7.82, 5.50, 4.88, 2.46];

%%%%% Plotting

createfigure_tput_mobihoc(Qth_plus_Contention_rho, Qth_plus_Contention_total_tput, ...
    CHAIN_1_rho, CHAIN_1_total_tput, ...
    QLB_CSMA_rho, QLB_CSMA_total_tput, ...
    DCF_rho, DCF_total_tput, case_code, filepath_tput);
createfigure_delay_mobihoc(Qth_plus_Contention_rho, Qth_plus_Contention_total_delay_s, ...
    CHAIN_1_rho, CHAIN_1_total_delay_s, ...
    QLB_CSMA_rho, QLB_CSMA_total_delay_s,...     
    DCF_rho, DCF_total_delay_s, case_code, filepath_delay);

createfigure_contention_percentage_mobihoc(Qth_plus_Contention_rho, Qth_plus_Contention_contention_time_percentage, ...
    CHAIN_1_rho, CHAIN_1_contention_time_percentage, ...
    QLB_CSMA_rho, QLB_CSMA_contention_time_percentage,...
    DCF_rho, DCF_contention_time_percentage, case_code, filepath_contention_percentage);

createfigure_collision_percentage_mobihoc(Qth_plus_Contention_rho, Qth_plus_Contention_collision_percentage, ...
    CHAIN_1_rho, CHAIN_1_collision_percentage, ...
    QLB_CSMA_rho, QLB_CSMA_collision_percentage,...    
    DCF_rho, DCF_collision_percentage, case_code, filepath_collision_percentage);


case '500-1-0-0-2groups-het-p=1'
%% Case: 500 nodes, 1 Chain, 0 DCF nodes, 2 groups, heterogeneous arrival rates
%% Qth-plus-Contention
Qth_plus_Contention_rho = [6.5, 6, 5.7, 5.5, 5, 4.5, 4, 3.5, 3];
Qth_plus_Contention_total_tput = [5770.38, 5544.4, 5318.18, 5141.27, 4730.57, 4260.26, 3785.95, 3312.68, 2839.19];
Qth_plus_Contention_total_delay_s = [7.954, 4.381, 3.547, 2.684, 1.044, 0.528, 0.591, 0.675, 0.786];
Qth_plus_Contention_contention_time_percentage = [5.16, 6.91, 7.89, 8.78, 13.95, 19.8, 20.28, 19.76 18.74];
Qth_plus_Contention_collision_percentage = [61.35, 59.48, 57.38, 55.41, 26.17, 9.81, 6.76, 4.99, 3.72];



%% Original CHAIN-1
CHAIN_1_rho = [6.5, 6, 5.5, 5, 4.5, 4, 3.7, 3.5, 3];
CHAIN_1_total_tput = [3983.73, 3974.91, 3957.48, 3998.44, 3919.39, 3733.43, 3510.86, 3319.01, 2851.85];
CHAIN_1_total_delay_s = [31.975, 27.359, 21.838, 14.120, 7.648, 2.037, 0.15, 0.00855, 0.00059];
CHAIN_1_contention_time_percentage = [36.26, 36.4, 36.68, 36.03, 37.29, 40.26, 43.79, 40.17, 28.12];
CHAIN_1_collision_percentage = [61.87, 61.34, 60.59, 59.96, 58.31, 54.65, 47.77, 22.53, 5.61];


%% QLB-CSMA
QLB_CSMA_rho = [6.5, 6, 5.5, 5, 4.5, 4, 3.6, 3.5, 3, 2.5];
QLB_CSMA_total_tput = [1206.95, 1210.06, 1204.56, 1207.24, 1212.09, 1229.53, 3416.78, 3328.22, 2846.49, 2373.42];
QLB_CSMA_total_delay_s = [72.542, 71.068, 69.456, 67.063, 64.377, 60.312, 5.36e-3, 4.82e-3, 3.29e-3, 2.59e-3];
QLB_CSMA_contention_time_percentage = [80.69, 80.64, 80.73, 80.68, 80.61, 75.44, 45.33, 46.75, 54.42, 61.78];
QLB_CSMA_collision_percentage = [75.93, 75.86, 75.97, 75.92, 75.82, 80.33, 6.35, 5.58, 3.34, 2.24];


%% DCF
DCF_rho = [6.5, 6, 5.5, 5, 4.5, 4, 3.5, 3.2, 3, 2.5];
DCF_total_tput = [1652.25, 1652.5, 1657.89, 1659.45, 1671.01, 1719.34, 2537.79, 3036.39, 2852.36, 2377.4];
DCF_total_delay_s = [65.845, 63.685, 61.365, 58.444, 54.114, 47.991, 21.031, 0.00114, 0.00062, 0.00039];
DCF_contention_time_percentage = [73.56, 73.56, 73.47, 73.45, 73.26, 72.49, 59.22, 31.15, 28.36, 23.39];
DCF_collision_percentage = [66.29, 66.29, 66.18, 66.14, 65.91, 64.9, 47.59, 8.69, 5.86, 2.86];

%%%%% Plotting

createfigure_tput_mobihoc(Qth_plus_Contention_rho, Qth_plus_Contention_total_tput, ...
    CHAIN_1_rho, CHAIN_1_total_tput, ...
    QLB_CSMA_rho, QLB_CSMA_total_tput, ...
    DCF_rho, DCF_total_tput, case_code, filepath_tput);
createfigure_delay_mobihoc(Qth_plus_Contention_rho, Qth_plus_Contention_total_delay_s, ...
    CHAIN_1_rho, CHAIN_1_total_delay_s, ...
    QLB_CSMA_rho, QLB_CSMA_total_delay_s,...     
    DCF_rho, DCF_total_delay_s, case_code, filepath_delay);

createfigure_contention_percentage_mobihoc(Qth_plus_Contention_rho, Qth_plus_Contention_contention_time_percentage, ...
    CHAIN_1_rho, CHAIN_1_contention_time_percentage, ...
    QLB_CSMA_rho, QLB_CSMA_contention_time_percentage,...
    DCF_rho, DCF_contention_time_percentage, case_code, filepath_contention_percentage);

createfigure_collision_percentage_mobihoc(Qth_plus_Contention_rho, Qth_plus_Contention_collision_percentage, ...
    CHAIN_1_rho, CHAIN_1_collision_percentage, ...
    QLB_CSMA_rho, QLB_CSMA_collision_percentage,...    
    DCF_rho, DCF_collision_percentage, case_code, filepath_collision_percentage);
end






