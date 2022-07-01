function plotMaxGrowth(mu_max,data,source_name,pos_color,pos_symbol)
% plotMaxGrowth
%   
%   Usage: plotMaxGrowth(mu_max,data,source_name,pos_color,pos_symbol)
%
%   Benjamin J. Sanchez, 2018-11-22
%

colors = [0    0    128             %dark blue
          0    0    255             %blue
          0    170  255             %light blue
          139  69   19              %saddlebrown
          0    128  0               %forest green
          0    255  0               %green
          128  255  0               %lime green
          128  128  0               %moss
          255  105  180             %pink
          255  255  0               %yellow
          191  0    255             %purple
          255  128  0               %orange
          0    0    0               %black
          255  0    0]./255;        %red
symbols = 'odssdddddddddddddddddds';

hold on
[M,N]       = size(mu_max);
data_resh   = reshape(data,M*N,1);
mu_max_resh = reshape(mu_max,M*N,1);
if nargin > 2
    %Fix color or symbol:
    if ~isempty(pos_color)
        colors = ones(size(colors)).*colors(pos_color,:);
    elseif ~isempty(pos_symbol)
        symbols = repmat(symbols(pos_symbol),[1 length(symbols)]);
    end
    %Sizes:
    s.la = 7;	%Axis label size
    s.nu = 7;	%Axis numbers size
    s.ma = 4;	%Marker size
    %Compute linear model:
    lmodel = fitlm(data_resh,mu_max_resh);
    R2     = lmodel.Rsquared.Ordinary;
    %Plot fit:
    axes = [0,0.3,0,0.5];
    plot(axes(1:2)',predict(lmodel,axes(1:2)'),'-k','LineWidth',1)
    text(axes(2) - 0.2,0.1,['R^{2} = ' num2str(round(R2,2))],'FontSize',s.la)
    %Add title:
    title(source_name)
else
    %Sizes:
    s.la = 15;	%Axis label size
    s.nu = 15;	%Axis numbers size
    s.ma = 8;	%Marker size
    %Compute mean error:
    mu_max_resh = mu_max_resh(data_resh > 0);
    data_resh   = data_resh(data_resh > 0);
    mean_error  = nanmean(abs(data_resh - mu_max_resh)./data_resh)*100;
    %Line y=x:
    axes = [0,0.5,0,0.5];
    plot(axes(1:2),axes(1:2),'-k','LineWidth',2)
    text(axes(2) - 0.3,0.1,['Mean error = ' num2str(round(mean_error,1)) '%'],'FontSize',s.la)
    %Axes:
    if max(max(data)) < 0.3
        xlabel('Experimental \mu_{max} [1/h]','FontSize',s.la)
    else
        xlabel('Rescaled experimental \mu_{max} [1/h]','FontSize',s.la)
    end
    ylabel('Predicted \mu_{max} [1/h]','FontSize',s.la)
end

%Plot:
for i = 1:M
    for j = 1:N
        plot(data(i,j),mu_max(i,j),symbols(i),'MarkerEdgeColor',colors(j,:), ...
             'MarkerFaceColor',colors(j,:),'MarkerSize',s.ma)
    end
end

%Other options:
axis square
xlim(axes(1:2))
ylim(axes(3:4))
set(gca,'FontSize',s.nu)
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gca,'XTick',axes(1):0.1:axes(2))
set(gca,'YTick',axes(3):0.1:axes(4))
box on
hold off

end
