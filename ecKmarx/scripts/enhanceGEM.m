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

%Provide your organism scientific name
org_name = 'kluyveromyces marxianus';
mkdir (['../models/' name])
%Convert model to RAVEN for easier visualization later on:
format short e
if strcmpi(toolbox,'COBRA')
    initCobraToolbox
    model = ravenCobraWrapper(model);
end

%Remove blocked rxns + correct model.rev:
cd change_model
[model,name,version] = preprocessModel(model,name,version);

%Retrieve kcats & MWs for each rxn in model:
cd ../get_enzyme_data
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,org_name);

%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);
[ecModel,modifications] = manualModifications(ecModel);

%Constrain model to batch conditions:
sigma    = 0.5;      
Ptot     = 0.3993;   %Assumed constant
gR_exp   = 0.44;     %[g/gDw h] Max batch gRate on minimal glucose media
c_source = 'D-glucose exchange (reversible)'; %Rxn name for the glucose uptake reaction
cd ../limit_proteins
[ecModel_batch,OptSigma] = getConstrainedModel(ecModel,c_source,sigma,Ptot,gR_exp,modifications,name);
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])

%Save output models:
cd ../../models
ecModel = saveECmodel(ecModel,toolbox,name,version);
ecModel_batch = saveECmodel(ecModel_batch,toolbox,[name '_batch'],version);
cd ../geckomat

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
