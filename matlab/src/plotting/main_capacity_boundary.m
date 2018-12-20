%% Scatter plots for the boundary of capacity region
clear;
clc;
filepath = '../figures/capacity-boundary/';
Xlimit_1 = [0 25];
Ylimit_1 = [0 25];
Xlimit_2 = [0 20];
Ylimit_2 = [0 20];

case_code = '500-1-0-0-2groups-p=1';

switch case_code
      
case '500-1-0-0-2groups-p=1'
Qth_plus_Qth_3_Contention_group1 = [23 21 19 17 15 13 11 10 9 8 6 4 2 1 0];
Qth_plus_Qth_3_Contention_group2 = [0 1 2 4 6 8 9 10 11 13 15 17 19 21 23];

Qth_plus_Qth_10_Contention_group1 = [24 22 20 18 16 15 14 13 11 10 9 8 7 5 3 1 0];
Qth_plus_Qth_10_Contention_group2 = [0 1 3 5 7 8 9 10 11 13 14 15 16 18 20 22 24];

Qth_plus_Qth_100_Contention_group1 = [24 23 22 20 18 16 14 13 12 11 10 9 7 5 4 2 1 0];
Qth_plus_Qth_100_Contention_group2 = [0 1 2 4 5 7 9 10 11 12 13 14 16 18 20 22 23 24];

CHAIN_group1 = [22 20 17 14 12 11 10 9 9 10 10 9 7 5 4 3 3 2 1 0 0];
CHAIN_group2 = [0 0 1 2 3 3 4 5 7 9 10 10 9 9 10 11 12 14 17 20 22];

QLB_CSMA_group1 = [14 12 10 9 8 7 6 5 4 3 1 0];
QLB_CSMA_group2 = [0 1 3 4 5 6 7 8 9 10 12 14];

DCF_group1 = [12 10 9 8 7 6 5 4 3 2 1 0];
DCF_group2 = [0 1 2 3 4 5 6 7 8 9 10 12];

%createfigure_capacity_boundary_QCHAIN(Qth_plus_Qth_3_Contention_group1, Qth_plus_Qth_3_Contention_group2,...
%                                      Qth_plus_Qth_10_Contention_group1, Qth_plus_Qth_10_Contention_group2,...
%                                      Qth_plus_Qth_100_Contention_group1, Qth_plus_Qth_100_Contention_group2);

createfigure_capacity_boundary_all(DCF_group1, DCF_group2,...
                                    QLB_CSMA_group1, QLB_CSMA_group2,...
                                      CHAIN_group1, CHAIN_group2,...
                                      Qth_plus_Qth_3_Contention_group1, Qth_plus_Qth_3_Contention_group2,...
                                      Qth_plus_Qth_10_Contention_group1, Qth_plus_Qth_10_Contention_group2,...
                                      Qth_plus_Qth_100_Contention_group1, Qth_plus_Qth_100_Contention_group2,...
                                      Xlimit_1, Ylimit_1);                                 

case '500-1-0-0-2groups-p=0.9'
Qth_plus_Qth_3_Contention_group1 = [14 13 11 10 9 8 6 5 4 3 2 1 0];
Qth_plus_Qth_3_Contention_group2 = [0 1 2 3 4 5 6 8 9 10 11 13 14];

Qth_plus_Qth_10_Contention_group1 = [17 16 14 12 11 10 9 8 7 6 5 3 1 0];
Qth_plus_Qth_10_Contention_group2 = [0 1 3 5 6 7 8 9 10 11 12 14 16 17];

Qth_plus_Qth_100_Contention_group1 = [18 17 15 13 12 11 10 9 8 7 6 5 3 1 0];
Qth_plus_Qth_100_Contention_group2 = [0 1 3 5 6 7 8 9 10 11 12 13 15 17 18];

CHAIN_group1 = [17 15 11 9 8 8 8 8 7 6 4 3 2 1 0];
CHAIN_group2 = [0 1 2 3 4 6 7 8 8 8 8 9 11 15 17];

QLB_CSMA_group1 = [12 11 10 9 8 7 6 5 4 3 2 1 0];
QLB_CSMA_group2 = [0 1 2 3 4 5 6 7 8 9 10 11 12];

DCF_group1 = [10 9 8 7 6 5 4 3 2 1 0];
DCF_group2 = [0 1 2 3 4 5 6 7 8 9 10];    
 

createfigure_capacity_boundary_all(DCF_group1, DCF_group2,...
                                      QLB_CSMA_group1, QLB_CSMA_group2,...
                                      CHAIN_group1, CHAIN_group2,...
                                      Qth_plus_Qth_3_Contention_group1, Qth_plus_Qth_3_Contention_group2,...
                                      Qth_plus_Qth_10_Contention_group1, Qth_plus_Qth_10_Contention_group2,...
                                      Qth_plus_Qth_100_Contention_group1, Qth_plus_Qth_100_Contention_group2,...
                                      Xlimit_2, Ylimit_2);   
end
%plot(Qth_plus_Qth_3_Contention_group1, Qth_plus_Qth_3_Contention_group2, '-bo');
%hold on
%plot(Qth_plus_Qth_10_Contention_group1, Qth_plus_Qth_10_Contention_group2, '-g+');
%plot(Qth_plus_Qth_100_Contention_group1, Qth_plus_Qth_100_Contention_group2, '-rs');
%{
Qth_plus_Qth_3_Contention_group1 = [23 22 21 20 19 18 17 16 15 15 14 13 13 12 11 10 9 8 8 7 6 6 5 4 4 3 2 1 1 0 0];
Qth_plus_Qth_3_Contention_group2 = [0 0 1 1 2 3 4 4 5 6 6 7 8 8 9 10 11 12 13 13 14 15 15 16 17 18 19 20 21 22 23];



Qth_plus_Qth_100_Contention_group1 = [24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 10 9 8 7 6 5 4 3 2 1 0];
Qth_plus_Qth_100_Contention_group2 = [0 1 2 3 4 5 6 7 8 9 10 10 11 12 13 14 15 16 17 18 19 20 21 22 24];



%}


%{
Qth_plus_Contention_Xmax = 580;
Qth_plus_Contention_y = int16([580 560 540 520 500 480 460 440 420 400 380 360 340 320 300 280 260 240 220 200 180 160 140 120 100 80 60 40 20 0]);
[Qth_plus_Contention_xout, Qth_plus_Contention_yout, Qth_plus_Contention_xout2, Qth_plus_Contention_yout2] = fill_interior(Qth_plus_Contention_Xmax, Qth_plus_Contention_y, Xlimit, Ylimit, Delta);

CHAIN_1_Xmax = 480;
CHAIN_1_y = int16([480 460 440 440 420 400 380 380 360 360 340 340 320 300 300 280 240 220 180 140 100 80 60 20 0]);
[CHAIN_1_xout, CHAIN_1_yout, CHAIN_1_xout2, CHAIN_1_yout2] = fill_interior(CHAIN_1_Xmax, CHAIN_1_y, Xlimit, Ylimit, Delta);

CHAIN_2_Xmax = 480;
CHAIN_2_y = int16([480 480 460 440 420 400 380 380 360 360 340 340 320 300 300 280 240 200 180 140 100 80 60 40 20]);
[CHAIN_2_xout, CHAIN_2_yout, CHAIN_2_xout2, CHAIN_2_yout2] = fill_interior(CHAIN_2_Xmax, CHAIN_2_y, Xlimit, Ylimit, Delta);

DCF_Xmax = 360;
DCF_y = int16([360 340 320 300 280 240 220 200 180 160 140 120 100 80 80 60 40 20 0]);
[DCF_xout, DCF_yout, DCF_xout2, DCF_yout2] = fill_interior(DCF_Xmax, DCF_y, Xlimit, Ylimit, Delta);


sz1 = 30;
sz2 = 15;
createfigure_cap(Qth_plus_Contention_xout, Qth_plus_Contention_yout, sz1, 'b', Qth_plus_Contention_xout2, Qth_plus_Contention_yout2, sz2, 'r', 'Autonomous CHAIN', filepath);
createfigure_cap(CHAIN_1_xout, CHAIN_1_yout, sz1, 'b', CHAIN_1_xout2, CHAIN_1_yout2, sz2, 'r', 'Original CHAIN-1', filepath);
createfigure_cap(CHAIN_2_xout, CHAIN_2_yout, sz1, 'b', CHAIN_2_xout2, CHAIN_2_yout2, sz2, 'r', 'Original CHAIN-2', filepath);
createfigure_cap(DCF_xout, DCF_yout, sz1, 'b', DCF_xout2, DCF_yout2, sz2, 'r', 'DCF', filepath);
%}
