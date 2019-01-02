%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotChemostats(protResults,condition)

%Exp data:
cd ./../exp_data
color  = sampleCVDmap(6);
ethanol = 'Ethanol production';
if strcmp(condition,'Temp')
    color    = color(6,:);
    cond     = 1;
    exp_data = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B2:G5');
    levels   = [30,33,36,38];
    x_label  = 'Temperature [°C]';
elseif strcmp(condition,'Osmo')
    color    = color(2,:);
    cond     = 2;
    levels   = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.3];
    x_label  = 'NaCl [M]';
    exp_data = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B7:G15');
elseif strcmp(condition,'EtOH')
    color    = color(4,:);
    cond     = 3;
    levels   = [0,20,40,60];
    x_label  = 'Etanol [g/L]';
    ethanol  = 'Ethanol consumption';
    exp_data = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B17:G20');
end
cd ./../results
exp_data = exp_data(1:length(levels),1:4);   %Glucose - O2 - CO2 - EtOH

%Model data:
mod_ori   = protResults.free{cond}(1:length(levels),2:5);
mod_wMC   = protResults.wMC{cond}(1:length(levels),2:5);
mod_wProt = protResults.wProt{cond}(1:length(levels),2:5);

%Exp data figure:
figure('position', [100,100,700,400])
plotData(levels,exp_data,'exp2',color,x_label,ethanol,'');

%Metabolic model:
figure('position', [100,100,700,400])
h = plotData(levels,exp_data,'exp',color,x_label,'','');
plotData(levels,mod_ori,'mod',color,x_label,ethanol,h);

%EC model - no proteomic data:
figure('position', [100,100,700,400])
h = plotData(levels,exp_data,'exp',color,x_label,'','');
plotData(levels,mod_wMC,'mod',color,x_label,ethanol,h);

%EC model - with proteomic data:
figure('position', [100,100,700,400])
h = plotData(levels,exp_data,'exp',color,'','');
plotData(levels,mod_wProt,'mod',color,x_label,ethanol,h);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = plotData(levels,data,type,color,x_label,ethanol,h_prev)

x    = levels;
y    = data;
ymax = 5*ceil(max(max(y))/5);

hold on
symbols = {'o','s','d','^'};
lines   = {'-','--',':','-.'};
for i = 1:length(data(1,:))
    if strcmp(type,'exp2')
        plot(x,y(:,i),[lines{i} 'o'],'Color',color,'MarkerSize',6, ...
               'MarkerFaceColor',color,'LineWidth',2);
    elseif strcmp(type,'exp')   
        h(i) = plot(x,y(:,i),symbols{i},'Color',color,'MarkerSize',6, ...
               'MarkerFaceColor',color,'LineWidth',2);
    else
        plot(x,y(:,i),lines{i},'Color',color,'LineWidth',2);
    end
end

%Other options:
text_size = 15;
xlim([x(1) x(end)])
set(gca,'XTick',x)
xlabel(x_label,'FontSize',text_size);
ylim(gca,[0 ymax])
ylabel('Flux [mmol/gDWh]','FontSize',text_size);
box on
set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
hold off

if ~isempty(ethanol)
    if isempty(h_prev)
        legend('Glucose uptake','O_{2} consumption','CO_{2} production', ...
               ethanol,'Location','northwest')
        legend('boxoff')
    else
        legend(h_prev,'Glucose uptake','O_{2} consumption','CO_{2} production', ...
               ethanol,'Location','northwest')
        legend('boxoff')
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
