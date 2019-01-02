%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function corrCond(stress,useP,names)

x = stress;
fit = zeros(length(names),2);
for i = 1:length(names)
    y   = useP(i,:);
    pos = ~isnan(y);
    if sum(pos) > 2
        lmodel   = fitlm(x(pos),y(pos));
        fit(i,1) = lmodel.Coefficients{2,1};
        fit(i,2) = sign(fit(i,1))*sqrt(lmodel.Rsquared.Ordinary);
    end
end

%Order from higher to lower and plot top 10:
[~,order] = sort(fit(:,1),'descend');
useP      = useP(order,:);
fit       = fit(order,:);
names     = names(order,:);
if length(names) > 100
    t     = 0.95;
    y_lab = 'm (slope of fit)';
else
    t     = 0.90;
    y_lab = 'Average slope of fit';
end

%Filter out if corr < t (in case of increase) of corr > -t (in case of
%decrease), and also if |slope| < 1e-3:
pos_inc = boolean((fit(:,2) >  t).*(abs(fit(:,1)) > 0.01));
pos_dec = boolean((fit(:,2) < -t).*(abs(fit(:,1)) > 0.01));

plotTopChanges(stress,fit,useP,'increasing',names,pos_inc,y_lab)
plotTopChanges(stress,fit,useP,'decreasing',names,pos_dec,y_lab)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotTopChanges(stress,fit,useP,direction,names,pos,y_lab)

if stress(4) == 38
    color = 'r';
    sname = 'temperature';
elseif stress(4) == 0.6
    color = 'b';
    sname = 'osmotic';
else
    color = [0 128 0]./255;
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
    pos_y = pos_y(1:n)+maxy/40;
    txt   = num2str(txt(1:n));
    txt   = strcat(txt,{'%'});
    text(pos_x,pos_y,txt,'HorizontalAlignment','center')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
