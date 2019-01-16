%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotStress(results,condition)

colors  = sampleCVDmap(6);
ethanol = '   Ethanol\newlineproduction';
if startsWith(condition,'REF')
    colors = [0,0,0];
    text   = '';
    cond   = [1,1];
elseif startsWith(condition,'Temp')
    colors = colors(6,:);
    text   = [' - ' condition(5:end) '°C'];
    if endsWith(condition,'33')
        cond = [1,2];
    elseif endsWith(condition,'36')
        cond = [1,3];
    elseif endsWith(condition,'38')
        cond = [1,4];
    end
elseif startsWith(condition,'Osmo')
    colors = colors(2,:);
    text  = [' - ' condition(5:end) ' M'];
    if endsWith(condition,'0.2')
        cond = [2,2];
    elseif endsWith(condition,'0.4')
        cond = [2,3];
    elseif endsWith(condition,'0.6')
        cond = [2,4];
    elseif endsWith(condition,'0.8')
        cond = [2,5];
    elseif endsWith(condition,'1.0')
        cond = [2,6];
    elseif endsWith(condition,'1.2')
        cond = [2,7];
    elseif endsWith(condition,'1.3')
        cond = [2,8];
    end
elseif startsWith(condition,'EtOH')
    colors  = colors(4,:);
    ethanol = '     Ethanol\newlineconsumption';
    text    = [' - ' condition(5:end) ' g/L'];
    if endsWith(condition,'20')
        cond = [3,2];
    elseif endsWith(condition,'40')
        cond = [3,3];
    elseif endsWith(condition,'60')
        cond = [3,4];
    end
end

%Exp data:
cd ./../exp_data
exp_data{1} = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B2:G5');     %Temp
exp_data{2} = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B7:G15');    %Osmotic Stress
exp_data{3} = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B17:G20');   %Ethanol
cd ./../results
exp_data = exp_data{cond(1)}(cond(2),1:4);   %Glucose - O2 - CO2 - EtOH

%Model predictions:
mod_free = results.free{cond(1)}(cond(2),2:5);
mod_wMC  = results.wMC{cond(1)}(cond(2),2:5);
if isfield(results,'wProt')
    mod_wProt = results.wProt{cond(1)}(cond(2),2:5);
    data      = [exp_data' mod_free' mod_wMC' mod_wProt'];
else
    data = [exp_data' mod_wMC' mod_free'];
end

%Plot data:
if isfield(results,'wProt')
    figure('position', [100,100,1200,600])
else
    figure('position', [100,100,800,500])
end
hold on
b = bar(data,'BarWidth',1);
legend(b,['Experimental data (0.1 1/h' text ')'],'Yeast8.0', ...
       'ecYeast8.0 - no proteomics', ...
       'ecYeast8.0 - with proteomics', ...
       'Location','northeast')

%Colors:
color(:,1) = linspace(1,colors(1),4)';
color(:,2) = linspace(1,colors(2),4)';
color(:,3) = linspace(1,colors(3),4)';
for i = 1:length(b)
    set(b(i),'FaceColor',color(i,:));
end

%Various options:
text_size = 15;
set(gca,'XTick',1:4,'XTickLabel',{'','','',''})
[hx,~] = format_ticks_v2(gca,{'    Glucose\newlineconsumption', 'O_{2} consumption', ...
                              'CO_{2} production', ethanol});
set(hx,'FontSize',text_size,'HorizontalAlignment','center')
ylims = get(gca,'YLim');
ylim([0 ylims(2)])
ylabel('Flux [mmol/gDWh]','FontSize',text_size);
legend('boxoff')
set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
box on
hold off

disp('Flux exchanges - metabolic model:')
checkExchange(results.mModel,results.fluxfree{cond(1)}(:,cond(2)))
disp('Flux exchanges - ecModel no Prot:')
checkExchange(results.ecModels_wMC{cond(1),cond(2)},results.fluxwMC{cond(1)}(:,cond(2)))
disp('Flux exchanges - ecModel wProt:')
checkExchange(results.ecModels_wProt{cond(1),cond(2)},results.fluxwProt{cond(1)}(:,cond(2)))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
