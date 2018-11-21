function [rangeGEM,rangeEC,indexes,stats] = FVA_ecModel(model,ecModel,tol,c_source,blockedMets)
% FVA_ecModel
%  
% This function goes through each of the rxns in a metabolic model and
% gets its flux variability range, then the rxn is mapped into an EC
% version of it to perform the correspondent variability analysis and
% finally compares and plots the cumulative flux variability distributions. 
%
%   model       MATLAB GEM structure
%   ecModel     MATLAB ecGEM structure
%   tol         numerical tolerance for a flux and variability rangeto be 
%               considered as zero
%   c_source    rxnName for the main carbon source uptake reaction (which
%               is fixed on the model according to the obtained value from 
%               a batch growth optimization with the ecModel).
%   blockedMets metNames of metabolites which secretion should be blocked
%   
%   rangeGEM    Distribution of variability ranges for the original GEM
%   rangeEC     Distribution of variability ranges for the original ecGEM
%   indexes     Indexes (in the original model) of the reactions for which 
%               a feasible variability range was obtained in both models
%   stats       Some statistics of the variability distributions
% 
% usage: [rangeDist,rangeDist_EC] = FVA_ecModel(model,ecModel,BlockFlag)
% 
% Ivan Domenzain.      Last edited: 2018-11-21

current       = pwd;
rangeGEM   = [];
indexes       = [];
rangeEC = [];
%Get the index for all the non-objective rxns in the original irrevModel
rxnsIndxs    = find(model.c~=1);
%Set minimal glucose media for ecModel
cd ../../scripts
[ecModel,pos] = changeMedia_batch(ecModel,[c_source ' (reversible)'],'Min');
%RAVEN built-in function
%irrevModel = convertToIrrev(model);
%Block glucose and oxygen production
cd (current)
if nargin>3
    irrevModel   = block_production(model,blockedMets,true);
    ecModel      = block_production(ecModel,blockedMets,true);
end
%Gets the optimal value for ecirrevModel and fixes the objective value to
%this for both irrevModels
[OptimalValue,ecFluxDist,ecModel]  = fixObjective(ecModel);
%Fix carbon source uptake rate for the ecModel on the original model
carbonUptake = ecFluxDist(pos(1));
disp([c_source ': ' num2str(carbonUptake)])
c_source           = find(strcmpi(model.rxnNames,c_source));
model.lb(c_source) = -carbonUptake;
[~,FluxDist,model] = fixObjective(model,OptimalValue);
%irrevModel.ub(c_source) = carbonUptake;
%irrevModel.lb(c_source) = 0.9999*carbonUptake;
%[~,FluxDist,irrevModel] = fixObjective(irrevModel,OptimalValue);

% Get the variability range for each of non-objective reactions in the
% original irrevModel
for i=1:length(rxnsIndxs) 
    indx        = rxnsIndxs(i);
    rxnID       = model.rxns(indx);
    %mappedIndxs = rxnMapping(rxnID,irrevModel,false);
    FixedValues = [];
    range       = MAXmin_Optimizer(model,indx,FixedValues,tol);
    %If max and min were feasible then the optimization proceeds with
    %the ecModel
    if ~isempty(range)
        if ~(range<tol & abs(FluxDist)<tol)
            mappedIndxs = rxnMapping(rxnID,ecModel,true);
            FixedValues = ecFluxDist(mappedIndxs);
            rangeEC     = MAXmin_Optimizer(ecModel,mappedIndxs,FixedValues,tol);
            if ~isempty(rangeEC)
                rangeGEM    = [rangeGEM; range];
                rangeEC  = [rangeEC; rangeEC];
                indexes        = [indexes; indx];
            end
        end
    end

    disp(['ready with #' num2str(i)])
end
%Plot FV cumulative distributions
distributions = {rangeGEM, rangeEC};
legends       = {'model', 'ecModel'};
titleStr      = 'Flux variability cumulative distribution';
[~, stats]    = plotCumDist(distributions,legends,titleStr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OptimalValue, optFluxDist, irrevModel] = fixObjective(irrevModel,priorValue)
% Optimize and fixes objective value for GEM
objIndx  = find(irrevModel.c~=0);
if nargin ==2
    irrevModel.lb(objIndx) = 0.9999*priorValue;
    irrevModel.ub(objIndx) = priorValue;
else
    sol = solveLP(irrevModel);
    irrevModel.lb(objIndx) = 0.9999*sol.x(objIndx);
    irrevModel.ub(objIndx) = sol.x(objIndx);
end
sol = solveLP(irrevModel);
if ~isempty(sol.f)
    OptimalValue = sol.x(objIndx);
    optFluxDist  = sol.x;
end
disp(['The optimal value is ' num2str(OptimalValue)])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_param, stats_param] = plotCumDist(E_parameter,legends,titlestr)
figure
for i=1:length(E_parameter)
    str{i} = horzcat(legends{i},' (',num2str(length(E_parameter{i})),...
        ' / ',num2str(median(E_parameter{i})), ')');
    [y_param(i), stats_param(i)] = cdfplot(E_parameter{i});
    title(titlestr)
    ylabel('Cumulative distribution','FontSize',30,'FontWeight','bold');
    xlabel('Variability range [mmol/gDw h]','FontSize',30,'FontWeight','bold');
    set(gca, 'XScale', 'log')
    hold on
end
legend(y_param,str);
end

