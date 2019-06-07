%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Csources_simulations
% 
% Ivan Domenzain. Last modified: 2019-06-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function Csources_simulations(ecModel_batch)
current = pwd;
colors = [0    0    128             %dark blue
          0    0    255             %blue
          0    170  255             %light blue
          139  69   19              %saddlebrown
          0    128  0               %forest green
          0    255  0               %green
          128  255  0               %lime green
          128  128  0               %moss
          255  105  180             %pink
          255  255  0               %yellow
          191  0    255             %purple
          255  128  0               %orange
          0    0    0               %black
          255  0    0]./255;        %red

figure
axis square
cd data
file_name  = 'Csources_growthRates.txt';
fID        = fopen(file_name);
data       = textscan(fID,'%s %f %f %s %s %f','delimiter','\t');
efe        = fclose('all'); 
count      = 1;
conditions = [];
gRates_exp = data{2};
stdDev     = data{3};
rxnNames   = data{4};
rxnIds     = data{5};
Ptot_ME    = data{6};
gIndex     = find(ecModel_batch.c);
gRates_sim = [];
error      = 0;
errors     = [];
legendStr  = cell(1,length(gRates_exp));
cd (current)
for i=1:length(gRates_exp)
    cd ../scripts
    c_source    = rxnNames{i};
    disp(c_source)
    c_source    = [c_source ' (reversible)'];
    [model,pos] = changeMedia_batch(ecModel_batch,c_source,'Min');
    protPos     = find(contains(model.rxnNames,'prot_'));
    %model.ub(protPos) = 0.6;
    solution = solveLP(model);
    %Sigma factor 
    sigma = 0.41;
    f     = 0.3827;
    if isempty(solution.x)
        gRate = 0;
        [components,shadowPrices] = findMissingComponents(model,1000);
        Ptot(i)                   = predictProteinCost(model,protPos,sigma,f);
    else
        disp(solution.f*-1)
        gRate = solution.x(gIndex);
        printFluxes(model,solution.x)
    end
    tempModel            = model;
    tempModel.lb(gIndex) = 0.9999*gRate;
    Ptot(i)              = predictProteinCost(tempModel,protPos,sigma,f);
    gRates_sim = [gRates_sim; gRate];
    error      = error + ((gRates_exp(i)-gRate)^2);
    errors(i)  = -((gRates_exp(i)-gRate)/gRates_exp(i));
    plot(gRates_exp(i),gRates_sim(i),'o','MarkerSize',8,'MarkerFaceColor',colors(i,1:3),'MarkerEdgeColor',colors(i,1:3))
    hold on
    legendStr{i} = strrep(c_source,' exchange (reversible)',''); 
end
error = sqrt(error/length(errors));
errorProt = (Ptot' - Ptot_ME)./Ptot_ME;
title('Max growth rate on different carbon sources','FontSize',30,'FontWeight','bold')
ylabel('\mu_{max} predicted [h^{-1}]','FontSize',30,'FontWeight','bold');
xlabel('\mu_{max} experimental [h^{-1}]','FontSize',30,'FontWeight','bold');
legend(legendStr)
%end
%--------------------------------------------------------------------------
function [components,shadowPrices] = findMissingComponents(model,perturbation)
%Function that searches missing components in a medium formulation for a
%model to become feasible
[~,uptkeIndxs] = getExchangeRxns(model);
newModel       = model; 
components     = {};
shadowPrices   = [];
for i=1:length(uptkeIndxs)
    index     = uptkeIndxs(i);
    tempModel = model;
    if tempModel.ub(index)== 0
        tempModel.ub(index) = perturbation;
        sol = solveLP(tempModel);
        if ~isempty(sol.x)
            disp(model.rxnNames{index})
            components   = [components; model.rxnNames(index)];
            shadowPrices = [shadowPrices; (-sol.f/perturbation)];
        end
    
    end
end
end
%--------------------------------------------------------------------------
function Ptot = predictProteinCost(model,protPos,sigma,f)
model.ub(protPos) = 1000;
model.c(:) = 0;
model.c(protPos(1:end-1)) = -1;
sol = solveLP(model);
if ~isempty(sol.x)
    Ptot = sol.x(protPos(end));
end
Ptot = Ptot/(sigma*f);
end