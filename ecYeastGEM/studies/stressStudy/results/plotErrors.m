%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotErrors(results)

%Calculate NGAM & errors:
Pcontent = zeros(14,3);
NGAM     = zeros(14,3);
errors   = zeros(14,3);
for i = 1:3
    for j = 1:9
        pos = 0;
        if ~isempty(results.ecModels_wProt{i,j})
            if i == 1 && j == 1
                pos = 1;
            elseif i == 1 && j == 2
                pos = 2;
            elseif i == 1 && j == 3
                pos = 3;
            elseif i == 1 && j == 4
                pos = 4;
            elseif i == 2 && j == 2
                pos = 5;
            elseif i == 2 && j == 3
                pos = 6;
            elseif i == 2 && j == 4
                pos = 7;
            elseif i == 2 && j == 5
                pos = 8;
            elseif i == 2 && j == 6
                pos = 9;
            elseif i == 2 && j == 7
                pos = 10;
            elseif i == 2 && j == 8
                pos = 11;
            elseif i == 3 && j == 2
                pos = 12;
            elseif i == 3 && j == 3
                pos = 13;
            elseif i == 3 && j == 4
                pos = 14;
            end
            
            if pos > 0
                Pcontent(pos,1) = results.Pmatched(i,j);
                Pcontent(pos,2) = results.Ppool(i,j);
                Pcontent(pos,3) = results.Ptot(i,j) - sum(Pcontent(pos,1:2));
                NGAM(pos,1)     = results.free{i}(j,end-1);
                NGAM(pos,2)     = results.wMC{i}(j,end-1);
                NGAM(pos,3)     = results.wProt{i}(j,end-1);
                errors(pos,1)   = results.free{i}(j,end);
                errors(pos,2)   = results.wMC{i}(j,end);
                errors(pos,3)   = results.wProt{i}(j,end);
            end
        end
    end
end

%Plot protein contents:
legend_text = {'Fraction in model, matched to proteomics', ...
               'Fraction in model, unmatched to proteomics', ...
               'Fraction not in model'};
plotBarPlot(Pcontent,'Protein content [g/gDW]',legend_text)

%Plot NGAM:
legend_text = {'Yeast8.0','ecYeast8.0 - no proteomic data', ...
               'ecYeast8.0 - with proteomic data'};
plotBarPlot(NGAM,'Fitted NGAM [mmol/gDWh]',legend_text)

%Spider plot for errors:
cd spider
figure('units', 'normalized', 'outerposition', [0 0.05 1 0.95]);
stress_names = {'  Reference',' 33�C','36�C','38�C','0.2 M','0.4 M','0.6 M ', ...
                '0.8 M ','1.0 M ','1.2 M','1.3 M','20 g/L','40 g/L',' 60 g/L'};
spider_plot(errors',stress_names',4,'LineStyle','-','LineWidth',2)
legend('Yeast8.0','ecYeast8.0 - no proteomics', ...
       'ecYeast8.0 - with proteomics','Location', 'southoutside');
legend('boxoff')
cd ..

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotBarPlot(data,y_label,legend_text)

figure('position', [50,50,1400,600])
colors = sampleCVDmap(6);
hold on

if contains(y_label,'content')
    b = bar(data,'stacked','BarWidth',0.5);
    set(b(1),'FaceColor',colors(1,:));
    set(b(2),'FaceColor',colors(3,:));
    set(b(3),'FaceColor',colors(5,:));
    max_e = 1;
    step  = 0.2;
else
    for i = 1:14
        b1 = bar(i-0.25,data(i,1),'BarWidth',0.25);
        b2 = bar(i,data(i,2),'BarWidth',0.25);
        b3 = bar(i+0.25,data(i,3),'BarWidth',0.25);
        if i == 1
            colors_i = [0,0,0];
        elseif i <= 4
            colors_i = colors(6,:);
        elseif i <= 11
            colors_i = colors(2,:);
        else
            colors_i = colors(4,:);
        end
        color(:,1) = linspace(1,colors_i(1),4)';
        color(:,2) = linspace(1,colors_i(2),4)';
        color(:,3) = linspace(1,colors_i(3),4)';
        set(b1,'FaceColor',color(2,:)) %Model
        set(b2,'FaceColor',color(3,:)) %ecModel_general
        set(b3,'FaceColor',color(4,:)) %ecModel_specific
    end
    b     = [b1 b2 b3];
    max_e = ceil(max(max(data))/10)*10;
    step  = 10;
end

%Plot separation lines:
plot([1.5,1.5],[0,max_e],'--k','LineWidth',1)
plot([4.5,4.5],[0,max_e],'--k','LineWidth',1)
plot([11.5,11.5],[0,max_e],'--k','LineWidth',1)

%Various options:
text_size = 15;
xlim([0.5,14.5])
ylim([0 max_e])
stress_names = {'Reference','33�C','36�C','38�C','0.2 M','0.4 M','0.6 M', ...
                '0.8 M','1.0 M','1.2 M','1.3 M','20 g/L','40 g/L','60 g/L'};
set(gca,'XTick',1:14,'XTickLabel',stress_names)
set(gca,'YTick',0:step:max_e)
ylabel(y_label,'FontSize',text_size);
legend(b,legend_text,'Location','north')
legend('boxoff')
set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
box on
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%