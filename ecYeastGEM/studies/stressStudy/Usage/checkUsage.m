%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prot,fluxes,singleCorr,loadings] = checkUsage(protResults)

delete(get(0,'Children'))   %deletes any present plot

REF_model = protResults.ecModel;
prot.use  = NaN(length(REF_model.enzymes),14);
prot.conc = NaN(length(REF_model.enzymes),14);
prot.useP = NaN(length(REF_model.enzymes),14);
not_prot  = cellfun(@isempty,strfind(REF_model.rxns,'prot'));
fluxes    = zeros(sum(not_prot),14);
cond      = 0;
for i = 1:length(protResults.ecModels(:,1))
    for j = 1:length(protResults.ecModels(1,:))
        ecModel = protResults.ecModels{i,j};
        if ~isempty(ecModel)
            cond = cond + 1;
            flux = protResults.fluxwProt{i}(:,j);
            fluxes(:,cond) = flux(not_prot);
            for k = 1:length(ecModel.enzymes)
                [usage,conc]      = getUsage(ecModel,ecModel.enzymes{k},flux);
                prot.use(k,cond)  = usage;
                prot.conc(k,cond) = conc;
                prot.useP(k,cond) = usage/conc*100;
            end
            %Check usage > 95%:
            u95 = sum(prot.useP(:,cond) > 95);
            disp(['condition #' num2str(cond) ': ' num2str(u95) ' enzymes with usage > 95%']);
            posEx = strcmp(ecModel.rxns,'prot_pool_exchange');
            if (ecModel.ub(posEx) - flux(posEx))/ecModel.ub(posEx) < 0.05
                disp('NOTE! unmeasured protein pool is limiting')
            end
        end
    end
end

%Construct table with enzyme properties: [kcat MW activity]
cd ..
kcats      = getKcats(REF_model);   %1/s
cd Usage
MWs        = REF_model.MWs;         %kDa = g/mmol
activs     = kcats./MWs*60;         %umol/mg/min
prot.props = [kcats MWs activs];

%Use gene names for any plot:
prot.codes = REF_model.enzymes;
prot.names = REF_model.geneNames;
for i = 1:length(prot.names)
    if strcmp(prot.names{i},'-')
        prot.names{i} = REF_model.enzymes{i};  %Change with UNIPROT code
    end
end

%Remove conditions with only NaNs in usage:
useP_means = nanmean(prot.useP,2);
only_NaNs  = isnan(useP_means);
prot.names = prot.names(~only_NaNs);
prot.codes = prot.codes(~only_NaNs);
prot.conc  = prot.conc(~only_NaNs,:);
prot.use   = prot.use(~only_NaNs,:);
prot.useP  = prot.useP(~only_NaNs,:);
prot.props = prot.props(~only_NaNs,:);

%Calculate averages and stds:
useP_means = nanmean(prot.useP,2);
useP_stds  = nanstd(prot.useP,0,2);

%Main metrics:
N = length(useP_means);
disp(['Average usage: ' num2str(mean(useP_means)) '%'])
disp(['Fraction of enzymes with usage < 10%: ' num2str(sum(useP_means<10)/N)])
disp(['Fraction of enzymes with usage > 50%: ' num2str(sum(useP_means>50)/N)])

%Histogram of usage:
figure('position',[0,0,600,600])
histPlot(useP_means,10,'Average usage [%]','enzymes')

%Plot top 10 used enzymes among all conditions:
yellow  = [254 210 24]./255;
topPlot(useP_means,useP_stds,prot.names,10,'Average usage [%]',[0 107],yellow,true)

%Correlation of percentual usage and 1.Kcats, 2.MWs, 3.Activities 4.Concentrations:
figure('position', [0,0,700,700])
subplot(2,2,1)
plotCorr(prot.props(:,1),useP_means,'Enzyme k_{cat} [1/s]',10.^(-2:2:6))
subplot(2,2,2)
plotCorr(prot.props(:,2),useP_means,'Enzyme MW [kDa]',0:50:250)
subplot(2,2,3)
plotCorr(prot.props(:,3),useP_means,'Enzyme s.a. [\mumol/(mg*min)]',10.^(-2:2:6))
subplot(2,2,4)
plotCorr(prot.conc,prot.useP,'Enzyme concentration [nmol/gDW]',10.^(-2:3))

%Correlation usage/condition:
Tlevels = [30 33 36 38];
Olevels = [0 0.2 0.4 0.6 0.8 1.0 1.2 1.3];
Elevels = [0 20 40 60];
figure('position', [0,0,1500,600])
subplot(1,3,1)
x_T = ones(size(prot.names))*Tlevels;
plotCorr(x_T,prot.useP(:,1:4),'Temperature stress [°C]',Tlevels)
subplot(1,3,2)
x_O = ones(size(prot.names))*Olevels;
plotCorr(x_O,prot.useP(:,[1 5:11]),'NaCl stress [M]',Olevels)
subplot(1,3,3)
x_E = ones(size(prot.names))*Elevels;
plotCorr(x_E,prot.useP(:,[1 12:14]),'Ethanol stress [g/L]',Elevels)
%2 small examples:
%figure('position', [0,0,400,400])
%pos = strcmp(prot.names,'ERG8');
%plotCorr(Tlevels,prot.useP(pos,1:4),'Temperature stress [°C]',Tlevels)
%figure('position', [0,0,400,400])
%pos = strcmp(prot.names,'TPS2');
%plotCorr(Tlevels,prot.useP(pos,1:4),'Temperature stress [°C]',Tlevels)

%Specific correlation usage/conditions:
corrCond(Tlevels,prot.useP(:,1:4),prot.names)
corrCond(Olevels,prot.useP(:,[1 5:11]),prot.names)
corrCond(Elevels,prot.useP(:,[1 12:14]),prot.names)

%Correlation usage/concentrations:
figure('position', [0,0,600,600])
plotCorr(prot.conc,prot.use,'Enzyme concentration [nmol/gDW]',10.^(-2:3))
%Small example:
%figure('position', [0,0,400,400])
%pos = strcmp(prot.names,'GCV2');
%plotCorr(prot.conc(pos,:),prot.use(pos,:),'Enzyme X',10.^(-2:3:5))

%Correlation usage/conc among conditions for each enzyme:
singleCorr = NaN(N,3);
for i = 1:N
    x = prot.conc(i,:);
    y = prot.use(i,:);
    x(isnan(x)) = 0;
    x = log10(x);
    y = log10(y);
    pos_x = ~isinf(x);
    pos_y = ~isinf(y);
    pos   = boolean(pos_x.*pos_y);
    if sum(pos_x) == 0 && sum(pos_y) > 0    %No usage in model of measured enzyme
        singleCorr(i,1) = 0;    %code
        singleCorr(i,2) = 0;    %slope
        singleCorr(i,3) = NaN;  %Pearson
    elseif sum(pos) > 2
        lmodel          = fitlm(x(pos),y(pos));
        singleCorr(i,1) = 1;
        singleCorr(i,2) = lmodel.Coefficients{2,1};
        singleCorr(i,3) = sign(singleCorr(i,2))*sqrt(lmodel.Rsquared.Ordinary);
    end
end
figure('position', [0,0,800,800])
subplot(2,2,1)
pos = ~isnan(singleCorr(:,3));
histPlot(singleCorr(pos,3),30,'r (Pearson correlation)','enzymes')
subplot(2,2,2)
pos = abs(singleCorr(:,2)) < 10;
histPlot(singleCorr(pos,2),30,'m (slope of fit)','enzymes')
subplot(2,2,3:4)
filter = boolean((singleCorr(:,3) > 0.95).*(singleCorr(:,2)<10));
topPlot(singleCorr(filter,2),[],prot.names(filter),10, ...
      'm (slope of fit)',[],yellow,false)
singleCorr(filter,2)

%PCA for conditions - logarithm:
loadings = PCAfigure(prot.conc,prot.useP,fluxes);

%Clustering + dendrograms + heat maps:
figure('position', [0,0,1500,800])
clusterStuff(log10(prot.conc),'conc','Enzymes',10,3,1,[])
clusterStuff(log10(prot.useP),'useP','Enzymes',10,3,2,[])
clusterStuff(log10(fluxes),'flux','Reactions',10,3,3,[])
figure('position', [0,0,800,800])
clusterStuff(log10(prot.useP),'useP','Enzymes',10,1,1,[])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [usage,conc] = getUsage(model,protein,flux)

pos = strcmp(model.rxns,['prot_' protein '_exchange']);
if sum(pos) == 0
    pos  = strcmp(model.rxns,['draw_prot_' protein]);
    conc = NaN;
else
    conc = model.ub(pos)*1e6;  %nmol/gDW
end
usage = flux(pos)*1e6;	%nmol/gDW

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loadings = PCAfigure(prot_concs,prot_useP,fluxes)

figure('position', [0,0,1500,500])
subplot(1,3,1)
[all_prot,loadings.prot] = PCAplot(prot_concs,'proteins');
subplot(1,3,2)
[all_useP,loadings.useP] = PCAplot(prot_useP,'usage');
subplot(1,3,3)
[all_flux,loadings.flux] = PCAplot(fluxes,'fluxes');

%Plot histogram with all data:
figure('position', [0,0,1500,500])
subplot(1,3,1);
histPlot(all_prot,[],'log10(concentration [nmol/gDW])','enzymes')
subplot(1,3,2);
histPlot(all_useP,[],'log10(usage [%])','enzymes')
subplot(1,3,3);
histPlot(all_flux,[],'log10(flux [mmol/gDWh])','reactions')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function histPlot(data,N,x_lab,y_lab)

if isempty(N)
    N = ceil(length(data)/20);
end

[counts,centers] = hist(data,N);
if ~isempty(strfind(x_lab,'slope'))
    bar(centers(counts > 0),counts(counts > 0),'c','BaseValue',0.7)
    set(gca,'yscale','log')
    y_lim   = [0.7 1e3];
    x_lim   = [-10,10];
    y_ticks = 10.^(0:3);
else
    bar(centers,counts,'c','BarWidth',1)
    y_lim   = [];
    x_lim   = [];
    y_ticks = [];
end
if N > 50
    set(get(gca,'child'),'EdgeColor','none');
end
setOptions(x_lab,x_lim,[],['Number of ' y_lab],y_lim,y_ticks)
if ~isempty(strfind(x_lab,'slope'))
    set(gca,'YTickLabel',{'1','10','100','1000'});
end
axis square

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotCorr(x,y,x_lab,x_tick)

%Replace any zero with NaN:
[m,n] = size(x);
for i = 1:n
    if isempty(strfind(x_lab,'stress'))
        pos = (x(:,i) == 0) + (y(:,i) == 0) > 0;
    else
        pos = (y(:,i) == 0) > 0;
    end
    x(pos,i) = NaN;
    y(pos,i) = NaN;
end

%Plot unfeasible region:
hold on
if max(max(y)) ~= 100 && m > 1
    patch([x_tick(1),x_tick(end),x_tick(1)], ...
          [x_tick(1),x_tick(end),x_tick(end)],[0.7 0.7 0.7]);
    text(x_tick(1)*2,x_tick(end-1),'Unfeasible region','FontSize',15)
end

%Plot data:
if n == 1
    colors = 'c';
elseif n == 14
    colors = getColors;
else
    if ~isempty(strfind(x_lab,'Temperature'))
        colors = [0 0 0;1 0 0;1 0 0;1 0 0];                         %red - temperature
    elseif ~isempty(strfind(x_lab,'NaCl'))
        colors = [0 0 0;0 0 1;0 0 1;0 0 1;0 0 1;0 0 1;0 0 1;0 0 1];	%blue - osmolarity
    else
        colors = [0 0 0;0 128 0;0 128 0;0 128 0]./255;              %forest green - ethanol
    end
end
x_r = [];
y_r = [];
if m == 1
    ms = 10;
    lw = 2;
else
    ms = 4;
    lw = 0.35;
end
if n == 14
    for i = [5:11 1 2 12 3 13 4 14]
        pos = isnan(x(:,i)+y(:,i)) == 0;
        plot(x(pos,i),y(pos,i),'o','MarkerEdgeColor',colors(i,:),'MarkerSize',ms,...
            'LineWidth',lw)
        x_r = [x_r;x(pos,i)];
        y_r = [y_r;y(pos,i)];
    end
else
    for i = 1:n
        pos = isnan(x(:,i)+y(:,i)) == 0;
        plot(x(pos,i),y(pos,i),'o','MarkerEdgeColor',colors(i,:),'MarkerSize',ms,...
            'LineWidth',lw)
        x_r = [x_r;x(pos,i)];
        y_r = [y_r;y(pos,i)];
    end
end

%Log scales and trendlines:
set(gca,'yscale','log')
if isempty(strfind(x_lab,'MW')) && isempty(strfind(x_lab,'stress'))
    set(gca,'xscale','log')
    lmodel = fitlm(log10(x_r),log10(y_r));
    plot(x_tick',10.^predict(lmodel,log10(x_tick')),'-m','LineWidth',lw)
else
    lmodel = fitlm(x_r,log10(y_r));
    plot(x_tick',10.^predict(lmodel,x_tick'),'-m','LineWidth',lw)
end
disp(['Fit for ' x_lab ': R^2 = ' num2str(lmodel.Rsquared.Ordinary)])

%Other options:
if max(max(y)) == 100
    y_tick = 10.^(-6:2:2);
    if n == 1
        y_lab  = 'Average usage [%]';
    else
        y_lab  = 'Enzyme usage [%]';
    end
elseif m == 1
    if ~isempty(strfind(x_lab,'stress'))
        if lmodel.Coefficients{2,1} > 0
            y_lab  = 'Enzyme X usage [%]';
            y_tick = ceil(10.^(1.4:0.05:1.55));
        else
            y_lab  = 'Enzyme Y usage [%]';
            y_tick = ceil(10.^(0.6:0.15:1.05));
        end
    else
        y_lab  = [];
        y_tick = 10.^(-2:3:5);
    end
else
    y_lab  = 'Enzyme usage [nmol/gDW]';
    y_tick = 10.^(-7:2:3);
end
setOptions(x_lab,[x_tick(1) x_tick(end)],x_tick, ...
            y_lab,[y_tick(1) y_tick(end)],y_tick)
axis square
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function colors = getColors

colors = [0    0    0               %black        - reference
          255  0    0               %red          - temperature
          255  0    0               %red          - temperature
          255  0    0               %red          - temperature
          0    0    255             %blue         - osmolarity
          0    0    255             %blue         - osmolarity
          0    0    255             %blue         - osmolarity
          0    0    255             %blue         - osmolarity
          0    0    255             %blue         - osmolarity
          0    0    255             %blue         - osmolarity
          0    0    255             %blue         - osmolarity
          0    128  0               %forest green - ethanol
          0    128  0               %forest green - ethanol
          0    128  0]./255;        %forest green - ethanol
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
