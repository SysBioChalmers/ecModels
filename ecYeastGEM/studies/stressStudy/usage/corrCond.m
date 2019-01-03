%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function corrProts = corrCond(var,usage,names)

[m,n] = size(var);
if m == 1 || n == 1
    percUsage = true;
else
    percUsage = false;
end

fit = NaN(length(names),2);
for i = 1:length(names)
    y = usage(i,:);
    if percUsage
        x = var;
    else
        x = var(i,:);
    end
    pos = ~isnan(y);
    % At least 3 data points for constructing linear fit:
    if sum(pos) >= 3
        lmodel   = fitlm(x(pos),y(pos));
        fit(i,1) = lmodel.Coefficients{2,1};
        fit(i,2) = lmodel.Rsquared.Ordinary;
    end
end

%Order from higher to lower slope:
[~,order] = sort(fit(:,1),'descend');
usage     = usage(order,:);
fit       = fit(order,:);
names     = names(order,:);

%Filters:
if length(names) > 100
    t     = 0.9;
    y_lab = 'Slope of fit';
else
    t     = 0.8;
    y_lab = 'Average slope of fit';
end
filter1 = fit(:,2) > t;             % R^2 > threshold
filter2 = ~isnan(sum(fit,2));       % No NaN values
filter3 = abs(fit(:,1)) > 1e-3;    % |slope| of at least > 1e-3
pos     = filter1.*filter2.*filter3 == 1;
if percUsage
    filter4 = max(usage,[],2) > 1;      % At least one value over 1% (for percentages)
    pos     = pos.*filter4 == 1;    
    %Distinguish between increasing usage and decreasing usage:
    pos_inc = pos.*(fit(:,1) > 0) == 1;
    pos_dec = pos.*(fit(:,1) < 0) == 1;
    %Plot top:
    plotTopChanges(var,fit,usage,'increasing',names,pos_inc,y_lab)
    plotTopChanges(var,fit,usage,'decreasing',names,pos_dec,y_lab)
    
else
    %"Volcano" plot:
    figure('position', [50,50,550,400])
    hold on
    hold all
    slope   = fit(filter2,1);
    R2      = fit(filter2,2);
    filter3 = filter3(filter2);
    plot(slope(filter3),R2(filter3),'o','MarkerEdgeColor','k', ...
         'MarkerFaceColor','r','MarkerSize',6,'LineWidth',0.5)
    plot(slope(~filter3),R2(~filter3),'o','MarkerEdgeColor','k', ...
         'MarkerFaceColor','k','MarkerSize',6,'LineWidth',0.5)
    limit_x = [-2,+2];
    limit_y = [t,t];
    plot(limit_x,limit_y,'--k','LineWidth',1)
    setOptions('Slope of fit',limit_x,[],'R^{2} of fit',[0 1],[])
    legend('Slope > 0.001','Slope < 0.001','Location','west')
    legend('boxoff')
    hold off
end

corrProts = names(pos);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotTopChanges(stress,fit,useP,direction,names,pos,y_lab)

colors = sampleCVDmap(6);
if stress(4) == 38
    color = colors(6,:);
    sname = 'temperature';
elseif stress(4) == 0.6
    color = colors(2,:);
    sname = 'osmotic';
else
    color = colors(4,:);
    sname = 'ethanol';
end

if strcmp(direction,'increasing')
    dir = 1;
else
    dir = -1;
end

if sum(pos) == 0
    disp(['No changes for ' direction ' ' sname ' stress!'])
end

maxy = max(dir*fit(pos,1));
if maxy > 1
    maxy = ceil(maxy);
elseif maxy > 0.1
    maxy = ceil(maxy*10)/10;
else
    maxy = ceil(maxy*100)/100;
end
topPlot(dir*fit(pos,1),[],names(pos),10,y_lab,[0 maxy],color,true)

if length(names) > 100
    %Indicate maximum usage in each case:
    txt   = max(useP,[],2);
    txt   = round(txt(pos)*10)/10;
    pos_y = dir*fit(pos,1);
    if strcmp(direction,'decreasing')
        txt   = flipud(txt);
        pos_y = flipud(pos_y);
    end
    
    %Show only selection:
    n     = min(10,sum(pos));
    pos_x = 1:n;
    pos_y = pos_y(1:n)+maxy/30;
    txt   = num2str(txt(1:n));
    txt   = strcat(txt,{'%'});
    text(pos_x,pos_y,txt,'HorizontalAlignment','center','FontSize',14)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
