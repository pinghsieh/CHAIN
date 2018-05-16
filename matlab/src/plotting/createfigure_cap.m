function createfigure_cap(X1, Y1, S1, C1, X2, Y2, S2, C2, myTitle, path)
%CREATEFIGURE2(X1, Y1, S1, C1, X2, Y2, S2, C2)
%  X1:  scatter x
%  Y1:  scatter y
%  S1:  scatter s
%  C1:  scatter c
%  X2:  scatter x
%  Y2:  scatter y
%  S2:  scatter s
%  C2:  scatter c

%  Auto-generated by MATLAB on 29-Jul-2017 21:03:29

% Create figure
figure1 = figure('Color',[1 1 1], 'position',[100,100,280,260]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create scatter
scatter(X1,Y1,S1,C1,'DisplayName','Achievable','MarkerFaceColor','flat',...
    'MarkerEdgeColor','none');

% Create scatter
scatter(X2,Y2,S2,C2,'DisplayName','Non-Achievable','Marker','x');

% Create xlabel
xlabel({'K_{Group 1} (packets/sec)'});

% Create ylabel
ylabel({'K_{Group 2} (packets/sec)'});

grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Helvetica Neue','FontSize',14,'LineWidth',1.5);

% Title
title(myTitle);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.461279461279461 0.772724802156409 0.442536381279377 0.128252788104089]);
savefig(strcat(path, myTitle, '.fig'));
saveas(gcf, strcat(path, myTitle, '.bmp'));

