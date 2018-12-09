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
proteins{4} = {'P37298','P33421'};      %Cytosolic one: Q00711 + P21801
protein_names{4} = 'Complex II';

% %NDI1/NDE1/NDE2:
% proteins{5} = {'P40215','Q07500','P32340'};
% protein_names{5} = 'NADH dehydrogenase';

%Settings depending on the condition:
if strfind(condition,'ref') > 0
    color  = [35  31  32        %Experimental
              115 115 119       %ecModel_general
              255 255 255]/255; %ecModel_specific
    text   = '';
    cond   = [1,1];
elseif strfind(condition,'temp') > 0
    color  = [237 31  36        %Experimental
              245 132 102       %ecModel_general
              255 255 255]/255; %ecModel_specific
    text   = [' - ' condition(5:end) '°C'];
    if strfind(text,'33') > 0
        cond   = [1,2];
    elseif strfind(text,'36') > 0
        cond   = [1,3];
    elseif strfind(text,'38') > 0
        cond   = [1,4];
    end
elseif strfind(condition,'osmo') > 0
    color  = [57  83  164       %Experimental
              128 134 193       %ecModel_general
              255 255 255]/255; %ecModel_specific
    text   = [' - ' condition(5:end) ' M'];
    if strfind(text,'0.2') > 0
        cond   = [2,2];
    elseif strfind(text,'0.4') > 0
        cond   = [2,3];
    elseif strfind(text,'0.6') > 0
        cond   = [2,4];
    elseif strfind(text,'0.8') > 0
        cond   = [2,5];
    elseif strfind(text,'1.0') > 0
        cond   = [2,6];
    elseif strfind(text,'1.2') > 0
        cond   = [2,7];
    elseif strfind(text,'1.3') > 0
        cond   = [2,8];
    end
elseif strfind(condition,'etoh') > 0
    color  = [13  129 64        %Experimental
              122 165 122       %ecModel_general
              255 255 255]/255; %ecModel_specific
    text    = [' - ' condition(5:end) ' g/L'];
    if strfind(text,'20') > 0
        cond   = [3,2];
    elseif strfind(text,'40') > 0
        cond   = [3,3];
    elseif strfind(text,'60') > 0
        cond   = [3,4];
    end
end

fluxeswMC   = results.fluxwMC{cond(1)}(:,cond(2));
fluxeswProt = results.fluxwProt{cond(1)}(:,cond(2));
values      = zeros(length(proteins),3);
for i = 1:length(proteins)
    for j = 1:length(proteins{i})
        %Usage in generic model:
        [usage,~]   = getUsage(results.ecModel,proteins{i}{j},fluxeswMC);
        values(i,2) = values(i,2) + usage;
        
        %Usage + concentration in specific model:
        [usage,ub]  = getUsage(results.ecModels{cond(1),cond(2)},proteins{i}{j},fluxeswProt);
        values(i,3) = values(i,3) + ub;
        values(i,1) = values(i,1) + usage;
    end
end
%Plot data:
figure('position', [50,50,700,300])
hold on
for i = 1:length(proteins)
    b1 = barh(i - 0.25,values(i,1),'BarWidth',0.25);
    b2 = barh(i,values(i,2),'BarWidth',0.25);
    b3 = barh(i + 0.25,values(i,3),'BarWidth',0.25);
    set(b1,'FaceColor',color(1,:)) %ecModel_specific
    set(b2,'FaceColor',color(2,:)) %ecModel_general
    set(b3,'FaceColor',color(3,:)) %Experimental
end

%Various options:
text_size = 15;
ylim([0.5,length(proteins) + 0.5])
set(gca,'YTick',1:length(proteins),'YTickLabel',protein_names)
xlabel('mg/gDW','FontSize',text_size);
legend([b3 b2 b1],['Protein concentration (0.1 1/h' text ')'], ...
                   'Protein usage - no proteomics', ...
                   'Protein usage - with proteomics', ...
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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
