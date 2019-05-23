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

%Constrain model to batch conditions:
cd ../limit_proteins
sigma    = 0.5;      
Ptot     = sumProtein(model);
gR_exp   = 0.278;     %[g/gDw h] Max batch gRate on minimal glucose media (CHASSY measurements)
c_source = 'D-glucose exchange (reversible)'; %Rxn name for the glucose uptake reaction

[ecModel_batch,OptSigma] = getConstrainedModel(ecModel,c_source,sigma,Ptot,gR_exp,modifications,name);
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])

%Save output models:
cd ../../models
ecModel = saveECmodel(ecModel,toolbox,name,version);
ecModel_batch = saveECmodel(ecModel_batch,toolbox,[name '_batch'],version);
cd ../geckomat

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
