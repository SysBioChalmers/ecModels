%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotProtein(results,condition)

%Complex V (ATP synthase):
proteins{1} = {'P07251','P00830','P38077','Q12165','P21306','P05626', ...
               'P30902','Q12349','P61829','P00854','P00856','Q06405', ...
               'P09457','P81450','P81451','Q12233','P81449'};
protein_names{1} = 'ATP synthase';

%Complex IV:
proteins{2} = {'P00401','Q01519','P32799','P00410','P00420','P04037', ...
               'P00427','P10174','P04039','P07255','P00425'};
protein_names{2} = 'Complex IV';

%Complex III:
proteins{3} = {'P00127','P00128','P00163','P07143','P07256','P07257', ...
               'P08067','P08525','P22289','P37299'};
protein_names{3} = 'Complex III';

%Complex II:
proteins{4} = {'P21801','P37298','P33421','P47052'};      %Cytosolic one: Q00711 + P21801
protein_names{4} = 'Complex II';

%Settings depending on the condition:
colors  = sampleCVDmap(6);
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
    text    = [' - ' condition(5:end) ' g/L'];
    if endsWith(condition,'20')
        cond = [3,2];
    elseif endsWith(condition,'40')
        cond = [3,3];
    elseif endsWith(condition,'60')
        cond = [3,4];
    end
end
color(:,1) = linspace(1,colors(1),4)';
color(:,2) = linspace(1,colors(2),4)';
color(:,3) = linspace(1,colors(3),4)';

fluxeswMC   = results.fluxwMC{cond(1)}(:,cond(2));
fluxeswProt = results.fluxwProt{cond(1)}(:,cond(2));
values      = zeros(length(proteins),3);
for i = 1:length(proteins)
    for j = 1:length(proteins{i})
        %Usage in generic model:
        [usage,~]   = getUsage(results.ecModels_wMC{cond(1),cond(2)},proteins{i}{j},fluxeswMC);
        values(i,2) = values(i,2) + usage;
        
        %Usage + concentration in specific model:
        [usage,ub]  = getUsage(results.ecModels_wProt{cond(1),cond(2)},proteins{i}{j},fluxeswProt);
        values(i,3) = values(i,3) + ub;
        values(i,1) = values(i,1) + usage;
    end
end
%Plot data:
figure('position', [50,50,600,500])
hold on
for i = 1:length(proteins)
    b1 = barh(i - 0.25,values(i,1),'BarWidth',0.25);
    b2 = barh(i,values(i,2),'BarWidth',0.25);
    b3 = barh(i + 0.25,values(i,3),'BarWidth',0.25);
    set(b1,'FaceColor',color(4,:)) %ecModel_specific
    set(b2,'FaceColor',color(3,:)) %ecModel_general
    set(b3,'FaceColor',color(1,:)) %Experimental
end

%Various options:
text_size = 15;
xlim([0,ceil(max(max(values))/5)*5])
ylim([0.5,length(proteins) + 0.5])
set(gca,'YTick',1:length(proteins),'YTickLabel',protein_names)
xlabel('mg/gDW','FontSize',text_size);
legend([b3 b2 b1],['Complex concentration (0.1 1/h' text ')'], ...
                   'Complex usage - no proteomics', ...
                   'Complex usage - with proteomics', ...
                   'Location','northeast')
legend('boxoff')
set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
box on
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [usage,ub] = getUsage(model,protein,flux)

MW  = model.MWs(strcmp(model.enzymes,protein));
pos = strcmp(model.rxns,['prot_' protein '_exchange']);
if sum(pos) == 0
    pos = strcmp(model.rxns,['draw_prot_' protein]);
end
usage = flux(pos)*MW*1000;
ub    = model.ub(pos)*MW*1000;
if isinf(ub)
    ub = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
