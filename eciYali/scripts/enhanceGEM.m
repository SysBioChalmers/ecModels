%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,ecModel_batch] = enhanceGEM(model,toolbox,name,version)
%
% Benjamin J. Sanchez & Ivan Domenzain. Last edited: 2019-02-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecModel,ecModel_batch] = enhanceGEM(model,toolbox,name,version)

if nargin < 3
    name    = '';
end
if nargin < 4
    version = '';
end
%Set lb for NGAM to 0.7 mmol ATP/gDw h (from S. cerevisiae)
model.lb(4) = 0.7;

%Provide your organism scientific name
org_name = 'yarrowia lipolytica';
mkdir (['../models/' name])

%Take out underscores from gene ids (better coverage in GECKO):
model.genes = strrep(model.genes,'YALI0_','YALI0');

% Remove compartment redundancy in metabolite ids/names (TODO: fix this in the RAVEN wrapper):
for i = 1:length(model.comps)
    comp = model.comps{i};
    model.mets = strrep(model.mets,['_' comp '[' comp ']'],['[' comp ']']);
    model.metNames = strrep(model.metNames,[' [' model.compNames{i} ']'],'');
end

%Convert model to RAVEN for easier visualization later on:
format short e
if isfield(model,'rules')
    initCobraToolbox
    model = ravenCobraWrapper(model);
end

%Remove blocked rxns + correct model.rev:
cd change_model
[model,name,version] = preprocessModel(model);

%Retrieve kcats & MWs for each rxn in model:
cd ../get_enzyme_data
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,org_name);

%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);
[ecModel,modifications] = manualModifications(ecModel);

%For a functional model, save upper bounds as +1000:
ecModel.ub(isinf(ecModel.ub)) = 1000;

%For ready to simulate model, modify bounds to minimal media:
cd ../kcat_sensitivity_analysis
[ecModel,~] = changeMedia_batch(ecModel,'D-glucose exchange (reversible)','Min',+1);

%Put back in underscores in gene ids (for visualization in caffeine):
ecModel = ravenCobraWrapper(ecModel);
ecModel.genes = strrep(ecModel.genes,'YALI0','YALI0_');
ecModel = ravenCobraWrapper(ecModel);

%Save output models:
cd ../../models
ecModel = saveECmodel(ecModel,toolbox,name,version);
cd ../geckomat

ecModel_batch = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
