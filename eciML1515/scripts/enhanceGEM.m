function [ecModel,ecModel_batch,version] = enhanceGEM(model,toolbox,name)
%enhanceGEM
%
%   Usage: [ecModel,ecModel_batch] = enhanceGEM(model,toolbox,name,version)
%
% Ivan Domenzain. Last edited: 2019-10-16
%

%Provide your organism scientific name
org_name = 'escherichia coli';
%Convert model to RAVEN for easier visualization later on:
format short e
if isfield(model,'rules')
    initCobraToolbox
    model = ravenCobraWrapper(model);
end
%Remove blocked rxns + correct model.rev:
cd change_model
[model,~,version] = preprocessModel(model,name);
%Retrieve kcats & MWs for each rxn in model:
cd ../get_enzyme_data
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,org_name);
%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);
[ecModel,modifications] = manualModifications(ecModel);
%Constrain model to batch conditions:
sigma    = 0.5;      %Optimized for glucose
Ptot     = 0.5793;      %Assumed constant
gR_exp   = 0.58;     %[g/gDw h] Max batch gRate on minimal glucose media
c_source = 'D-glucose exchange (reversible)'; %Rxn name for the glucose uptake reaction
cd ../limit_proteins
mkdir (['../../models/' name]);
[ecModel_batch,OptSigma] = getConstrainedModel(ecModel,c_source,sigma,Ptot,gR_exp,modifications,name);
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])
%Save output models:
cd ../../models
ecModel = saveECmodel(ecModel,toolbox,name,version);
ecModel_batch = saveECmodel(ecModel_batch,toolbox,[name '_batch'],version);
cd ../geckomat
end