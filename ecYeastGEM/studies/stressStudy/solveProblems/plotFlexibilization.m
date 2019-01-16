%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotFlexibilization(protResults,list,option)

list{2,1} = [];
list{3,1} = [];
ecModel   = protResults.ecModels_wProt{1,1};
data      = zeros(length(ecModel.enzymes),14);
n         = 0;
for i = 1:length(list(1,:))
    for j = 1:length(list(:,1))
        if ~isempty(list{j,i})
            n = n + 1;
            for k = 1:length(list{j,i}(:,1))
                pos = strcmp(ecModel.enzNames,list{j,i}{k,1});
                if list{j,i}{k,3} > list{j,i}{k,2}
                    switch option
                        case 'FC'
                            data(pos,n) = log10(list{j,i}{k,3}/list{j,i}{k,2});
                        case 'mass'
                            diff_mass   = (list{j,i}{k,3} - list{j,i}{k,2})*ecModel.MWs(pos);
                            data(pos,n) = diff_mass/protResults.Ptot(j,i)*100;
                    end
                end
            end
        end
    end
end

%Filter out enzymes with no flexibilization:
enzymes = ecModel.enzNames;
enzymes(sum(data,2) == 0) = [];
data(sum(data,2) == 0,:)  = [];

%Sort enzymes by order of appearance:
is_flex = data > 0;
order   = zeros(size(enzymes));
k       = 1;
for i = 1:length(data)
    pos = find(is_flex(:,i));
    for j = 1:length(pos)
        order(k) = pos(j);
        is_flex(pos(j),:) = 0;
        k = k + 1;
    end
end
enzymes = enzymes(order);
data    = data(order,:);

%Plot data:
figure('position',[50,50,1400,600])
hold on
text_size = 15;
cmap = parula(length(enzymes));
switch option
    case 'FC'
        b = bar(data',1.2);
        ylabel('log10(increase)','FontSize',text_size);
        max_y = ceil(max(max(data)));
    case 'mass'
        b = bar(data','stacked');
        ylabel('Added mass [%]','FontSize',text_size);
        disp(['max added mass: ' num2str(max(sum(data)))])
        max_y = ceil(max(sum(data))/5)*5;
end
for i = 1:length(enzymes)
    set(b(i),'FaceColor',cmap(i,:));
end
legend(b,enzymes','Location','NorthEast')
legend('boxoff')

%Plot separation lines:
plot([1.5,1.5],[0,max_y],'--k','LineWidth',1)
plot([4.5,4.5],[0,max_y],'--k','LineWidth',1)
plot([11.5,11.5],[0,max_y],'--k','LineWidth',1)

%Other options:
xlim([0.5,14.5])
ylim([0,max_y])
stress_names = {'Reference','33°C','36°C','38°C','0.2 M','0.4 M','0.6 M', ...
                '0.8 M','1.0 M','1.2 M','1.3 M','20 g/L','40 g/L','60 g/L'};
set(gca,'XTick',1:14,'XTickLabel',stress_names)
set(gca,'YTick',1:max_y)
xlabel('Dilution rate [1/h]','FontSize',text_size);
box on
set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
