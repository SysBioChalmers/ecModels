%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotStress(results,condition)

%Exp data:
cd ../exp_data
ethanol = '   Ethanol\newlineproduction';
no_std  = false;
if strfind(condition,'ref') > 0
    color  = [255 255 255       %Experimental
              210 209 212       %Model
              115 115 119       %ecModel_general
              35  31  32 ]/255; %ecModel_specific
    text   = '';
    values = 'K3:R6';     %Ref
    cond   = [1,1];
elseif strfind(condition,'temp') > 0
    color  = [255 255 255       %Experimental
              252 211 193       %Model
              245 132 102       %ecModel_general
              237 31  36 ]/255; %ecModel_specific
    text   = [' - ' condition(5:end) '°C'];
    if strfind(text,'33') > 0
        values = 'K7:R11';     %Temp33deg
        cond   = [1,2];
    elseif strfind(text,'36') > 0
        values = 'K12:R16';     %Temp36deg
        cond   = [1,3];
    elseif strfind(text,'38') > 0
        values = 'K17:R19';     %Temp38deg
        cond   = [1,4];
    end
elseif strfind(condition,'osmo') > 0
    color  = [255 255 255       %Experimental
              205 206 232       %Model
              128 134 193       %ecModel_general
              57  83  164]/255; %ecModel_specific
    text   = [' - ' condition(5:end) ' M'];
    if strfind(text,'0.2') > 0
        values = 'K20:R22';     %Osm0.2M
        cond   = [2,2];
    elseif strfind(text,'0.4') > 0
        values = 'K23:R25';     %Osm0.4M
        cond   = [2,3];
    elseif strfind(text,'0.6') > 0
        values = 'K26:R28';     %Osm0.6M
        cond   = [2,4];
    elseif strfind(text,'0.8') > 0
        values = 'K29:R31';     %Osm0.8M
        cond   = [2,5];
    elseif strfind(text,'1.0') > 0
        values = 'B12:G12';     %Osm1.0M
        cond   = [2,6];
        no_std = true;
    elseif strfind(text,'1.2') > 0
        values = 'B13:G13';     %Osm1.2M
        cond   = [2,7];
        no_std = true;
    elseif strfind(text,'1.3') > 0
        values = 'B14:G14';     %Osm1.3M
        cond   = [2,8];
        no_std = true;
    end
elseif strfind(condition,'etoh') > 0
    color  = [255 255 255       %Experimental
              204 216 200       %Model
              122 165 122       %ecModel_general
              13  129 64 ]/255; %ecModel_specific
    ethanol = '     Ethanol\newlineconsumption';
    text    = [' - ' condition(5:end) ' g/L'];
    if strfind(text,'20') > 0
        values = 'K32:R34';     %EtOH20g/L
        cond   = [3,2];
    elseif strfind(text,'40') > 0
        values = 'K35:R37';     %EtOH40g/L
        cond   = [3,3];
    elseif strfind(text,'60') > 0
        values = 'K38:R40';     %EtOH60g/L
        cond   = [3,4];
    end
end
if no_std
    exp_data = xlsread('Sce_stress_chemostats_merged.xlsx',1,values);
    exp_mean = exp_data(:,1:4);   %Glucose - O2 - CO2 - EtOH
else
    exp_data = xlsread('Sce_stress_chemostats.xlsx','Raw_data',values);
    exp_mean = abs(mean(exp_data(:,[1 7 8 4]),1));   %Glucose - O2 - CO2 - EtOH
end
cd ../results

%Model predictions:
mod_free = results.free{cond(1)}(cond(2),2:5);
mod_wMC  = results.wMC{cond(1)}(cond(2),2:5);
if isfield(results,'wProt')
    mod_wProt = results.wProt{cond(1)}(cond(2),2:5);
    data      = [exp_mean' mod_free' mod_wMC' mod_wProt'];
else
    data = [exp_mean' mod_wMC' mod_free'];
end

%Plot data:
if isfield(results,'wProt')
    figure('position', [100,100,1200,600])
else
    figure('position', [100,100,800,500])
end
hold on
b = bar(data,'BarWidth',1);
if isfield(results,'wProt')
    if strfind(condition,'etoh') > 0
        legend(b,['Experimental data (0.1 1/h' text ')'],'Yeast7', ...
                 'ecYeast7 - no proteomics', ...
                 'ecYeast7 - with proteomics', ...
                 'Location','northeast')
    else
        legend(b,['Experimental data (0.1 1/h' text ')'],'Yeast7', ...
                 'ecYeast7 - no proteomics', ...
                 'ecYeast7 - with proteomics', ...
                 'Location','northwest')
    end
else
    color = [1  1        0      %Yellow
             1  128/255  0      %Orange
             1  0        0];    %Red
    legend(b,['Experimental data (0.1 1/h' text ')'], ...
             'ecYeast7 model', ...
             'Yeast7 model', ...
             'Location','northwest')
end

for i = 1:length(b)
    set(b(i),'FaceColor',color(i,:));
end

%Plot stds for experimental data:
if ~no_std
    exp_std = std(exp_data(:,[1 7 8 4]),1);         %Glucose - O2 - CO2 - EtOH
    pos     = get(b,'XData');
    ofs     = get(b,'XOffset');
    errorbar(pos{1}+ofs{1},exp_mean,exp_std,'k','linestyle','none')
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
disp('Flux exchanges - ecYeast7 no Prot:')
checkExchange(results.ecModels{cond(1),cond(2)},results.fluxwMC{cond(1)}(:,cond(2)))
disp('Flux exchanges - ecYeast7 wProt:')
checkExchange(results.ecModels{cond(1),cond(2)},results.fluxwProt{cond(1)}(:,cond(2)))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
