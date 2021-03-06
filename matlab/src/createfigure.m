function createfigure(X1, YMatrix1)
%CREATEFIGURE1(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 28-Feb-2017 10:25:33

% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',2,'Parent',axes1);
s = size(YMatrix1, 1);
if s >= 1
set(plot1(1),'DisplayName','Node 1','Color',[0 0 1]);
end
if s >= 2
set(plot1(2),'DisplayName','Node 2','Color',[1 0 1]);
end
if s >= 3
set(plot1(3),'DisplayName','Node 3',...
    'Color',[0.980392156862745 0.588235294117647 0.196078431372549]);
end
if s >= 4
set(plot1(4),'DisplayName','Node 4','Color',[1 0 0]);
end
if s >= 5
set(plot1(5),'DisplayName','Node 5','Color',[0 0.784313725490196 0]);
end

% Create xlabel
xlabel({'Time'});

% Create ylabel
ylabel({'Queue Length'});

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',18,'LineWidth',2,'XGrid','on','YGrid','on');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.156392585598693 0.645238095238096 0.166964285714286 0.254529389880954]);

