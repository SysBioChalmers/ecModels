%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function topPlot(means,stds,names,N,y_lab,y_lim,color,new)

if isempty(means)
    return
end

%Show only one ATPsynt value:
ATP_pos = ~cellfun(@isempty,strfind(names,'ATP')) + strcmp(names,'OLI1') + strcmp(names,'TIM11') == 1;
if sum(ATP_pos) > 1
    names{find(ATP_pos,1)}   = 'ATPsynt';
    ATP_pos(find(ATP_pos,1)) = false;
    names = names(~ATP_pos);
    means = means(~ATP_pos);
    if ~isempty(stds)
        stds = stds(~ATP_pos);
    end
end

%Replace NaN with zeros:
means(isnan(means)) = 0;

%Order from higher sum to lower:
[~,order] = sort(means,'descend');
means     = means(order);
names     = names(order);

maxx = min(N,length(names));

%Plot top 10:
if new
    figure('position', [0,0,1400*(maxx+5)/20,600])
end
hold on
bar(means(1:maxx),'FaceColor',color,'BarWidth',0.7);

%Standard deviation:
if ~isempty(stds)
    stds(isnan(stds)) = 0;
    stds              = stds(order);
    errorbar(1:maxx,means(1:maxx),stds(1:maxx),'k','linestyle','none')
end

setOptions([],[0 maxx+1],1:min(N,length(names)),y_lab,y_lim,[])
set(gca,'XTickLabel',names(1:maxx),'FontSize',15)
if length(names{1}) > 6
    xticklabel_rotate([],45)
end
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
