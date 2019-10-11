%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prot,fluxes] = checkUsage(protResults)

delete(get(0,'Children'))   %deletes any present plot

REF_model = protResults.ecModels_wMC{1,1};
prot.use  = NaN(length(REF_model.enzymes),14);
prot.conc = NaN(length(REF_model.enzymes),14);
prot.useP = NaN(length(REF_model.enzymes),14);
not_prot  = ~contains(REF_model.rxns,'prot');
fluxes    = zeros(sum(not_prot),14);
cond      = 0;
protResults.ecModels_wProt{2,1} = '';
protResults.ecModels_wProt{3,1} = '';
for i = 1:length(protResults.ecModels_wProt(:,1))
    for j = 1:length(protResults.ecModels_wProt(1,:))
        ecModel = protResults.ecModels_wProt{i,j};
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
cd usage
MWs        = REF_model.MWs;         %kDa = g/mmol
activs     = kcats./MWs*60;         %umol/mg/min
prot.props = [kcats MWs activs];

%Use gene names for any plot:
prot.codes = REF_model.enzymes;
prot.names = REF_model.enzNames;
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
topPlot(useP_means,useP_stds,prot.names,22,'Average usage [%]',[0 107],yellow,true)

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
plotCorr(x_O,prot.useP(:,[1,5:11]),'NaCl stress [M]',Olevels)
subplot(1,3,3)
x_E = ones(size(prot.names))*Elevels;
plotCorr(x_E,prot.useP(:,[1,12:14]),'Ethanol stress [g/L]',Elevels)

%Specific correlation usage/conditions:
corrCond(Tlevels,prot.useP(:,1:4),prot.names);
corrCond(Olevels,prot.useP(:,[1 5:11]),prot.names);
corrCond(Elevels,prot.useP(:,[1 12:14]),prot.names);

%Correlation usage/concentrations:
figure('position', [0,0,600,600])
plotCorr(prot.conc,prot.use,'Enzyme concentration [nmol/gDW]',10.^(-2:3))
%TODO: 10^6 fmol/gDW = 1 nmol/gDW
prot.corrProts = corrCond(prot.conc,prot.use,prot.names);

%PCA for conditions:
PCAfigure(prot.conc,prot.use,fluxes);

%Clustering + dendrograms + heat maps:
figure('position', [0,0,1500,800])
clusterStuff(log10(prot.conc),'conc','Enzymes',10,3,1,[])
clusterStuff(log10(prot.useP),'useP','Enzymes',10,3,2,[])
clusterStuff(log10(fluxes),'flux','Reactions',10,3,3,[])
figure('position', [0,0,800,800])
clusterStuff(log10(prot.useP),'useP','Enzymes',10,1,1,[])

%Tryptophan story:
TRP_names = {'TRP1','TRP2','TRP3','TRP4','TRP5'};
TRP_table = zeros(length(TRP_names),4);
for i = 1:length(TRP_names)
    pos = strcmp(prot.names,TRP_names{i});
    if sum(pos) == 1
        TRP_table(i,:) = prot.useP(pos,1:4);
    end
end
map = flipud(sampleCVDmap(100));
colormap(map)
imagesc(TRP_table)
colorbar

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
