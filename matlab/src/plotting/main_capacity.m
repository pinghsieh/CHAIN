%% Scatter plots for capacity region
clear;
clc;
filepath = '../figures/capacity/';
Xlimit = 600;
Ylimit = 600;
Delta = 20;

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
