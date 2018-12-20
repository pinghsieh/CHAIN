function createfigure_capacity_boundary_all(X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5, X6, Y6, my_xlimit, my_ylimit)
%CREATEFIGURE2(X1, Y1, X2, Y2, X3, Y3)
%  X1:  vector of x data
%  Y1:  vector of y data
%  X2:  vector of x data
%  Y2:  vector of y data
%  X3:  vector of x data
%  Y3:  vector of y data

%  Auto-generated by MATLAB on 17-Dec-2018 16:59:26

% Create figure
figure1 = figure('Color',[1 1 1],'position',[100,100,450,450]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
plot(X1,Y1,'DisplayName','DCF','MarkerSize',10,'Marker','x',...
    'LineWidth',2.5,...
    'LineStyle',':',...
    'Color',[1 0 0]);

% Create plot
plot(X2,Y2,'DisplayName','QLB-CSMA','MarkerSize',10,'Marker','v',...
    'LineWidth',2.5,...
    'Color',[0 1 0]);

% Create plot
plot(X3,Y3,'DisplayName','CHAIN','MarkerSize',10,'Marker','^',...
    'LineWidth',2.5,...
    'Color',[0 0 1]);

% Create plot
plot(X4,Y4,'DisplayName','Q-CHAIN, Q*=3','MarkerSize',10,...
    'Marker','o',...
    'LineWidth',2.5,...
    'LineStyle','--',...
    'Color',[1 0.7 0]);

% Create plot
plot(X5,Y5,'DisplayName','Q-CHAIN, Q*=10','MarkerSize',11,...
    'Marker','square',...
    'LineWidth',2.5,...
    'LineStyle','--',...
    'Color',[1 0 1]);

% Create plot
plot(X6,Y6,'DisplayName','Q-CHAIN, Q*=100','MarkerSize',10,...
    'Marker','+',...
    'LineWidth',2.5,...
    'LineStyle','--',...
    'Color',[0 1 1]);

% Create ylabel
ylabel({'K_{Group2} (packets/sec)'});

% Create xlabel
xlabel({'K_{Group1}(packets/sec)'});

box(axes1,'on');
grid(axes1,'on');
xlim(my_xlimit);
ylim(my_ylimit);
% Set the remaining axes properties
set(axes1,'FontSize',20,'LineWidth',1.5);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.533035714285714 0.716666666666667 0.34375 0.182142857142859]);
