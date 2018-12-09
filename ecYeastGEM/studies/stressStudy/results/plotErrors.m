%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotErrors(results)

%Calculate errors:
errors  = zeros(14,3);
for i = 1:3
    for j = 1:9
        pos = 0;
        if ~isempty(results.ecModels{i,j})
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
            
            %Retrieve error:
            if pos > 0
                errors(pos,1) = results.free{i}(j,end);
                errors(pos,2) = results.wMC{i}(j,end);
                errors(pos,3) = results.wProt{i}(j,end);
            end
        end
    end
end

%Plot data:
figure('position', [50,50,1400,600])
hold on
for i = 1:14
    b1 = bar(i-0.25,errors(i,1),'BarWidth',0.25);
    b2 = bar(i,errors(i,2),'BarWidth',0.25);
    b3 = bar(i+0.25,errors(i,3),'BarWidth',0.25);
    if i == 1
        set(b1,'FaceColor',[210 209 212]/255) %Model
        set(b2,'FaceColor',[115 115 119]/255) %ecModel_general
        set(b3,'FaceColor',[35  31  32 ]/255) %ecModel_specific
    elseif i <= 4
        set(b1,'FaceColor',[252 211 193]/255) %Model
        set(b2,'FaceColor',[245 132 102]/255) %ecModel_general
        set(b3,'FaceColor',[237 31  36 ]/255) %ecModel_specific
    elseif i <= 11
        set(b1,'FaceColor',[205 206 232]/255) %Model
        set(b2,'FaceColor',[128 134 193]/255) %ecModel_general
        set(b3,'FaceColor',[57  83  164]/255) %ecModel_specific
    else
        set(b1,'FaceColor',[204 216 200]/255) %Model
        set(b2,'FaceColor',[122 165 122]/255) %ecModel_general
        set(b3,'FaceColor',[13  129 64 ]/255) %ecModel_specific
    end
    
end

%Plot separation lines:
max_e = ceil(max(max(errors))/10)*10;
plot([1.5,1.5],[0,max_e],'--k','LineWidth',1)
plot([4.5,4.5],[0,max_e],'--k','LineWidth',1)
plot([11.5,11.5],[0,max_e],'--k','LineWidth',1)

%Various options:
text_size = 15;
xlim([0.5,14.5])
ylim([0 max_e])
stress_names = {'Reference','33°C','36°C','38°C','0.2 M','0.4 M','0.6 M', ...
                '0.8 M','1.0 M','1.2 M','1.3 M','20 g/L','40 g/L','60 g/L'};
set(gca,'XTick',1:14,'XTickLabel',stress_names)
set(gca,'YTick',0:10:max_e)
ylabel('Average error in fitting [%]','FontSize',text_size);
legend([b1 b2 b3],'Yeast7','ecYeast7 - no proteomic data', ...
       'ecYeast7 - with proteomic data','Location','north')
legend('boxoff')
set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
box on
hold off

%Spider plot:
cd spider
figure('units', 'normalized', 'outerposition', [0 0.05 1 0.95]);
stress_names = {'  Reference',' 33°C','36°C','38°C','0.2 M','0.4 M','0.6 M ', ...
                '0.8 M ','1.0 M ','1.2 M','1.3 M','20 g/L','40 g/L',' 60 g/L'};
spider_plot(log(errors'),stress_names',4,'LineStyle','-','LineWidth',2)
legend('Yeast7','ecYeast7 - no proteomic data', ...
       'ecYeast7 - with proteomic data','Location', 'southoutside');
legend('boxoff')
cd ..

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
