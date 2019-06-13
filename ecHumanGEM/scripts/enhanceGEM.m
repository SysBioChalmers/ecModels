function [ecModel,ecModel_batch] = enhanceGEM(model,cellName)
%enhanceGEM_cellLine
% 
% Function that loads a cell-line or tissue specific human metabolism model
% and enhances it with enzyme-constraints with the use of the GECKO
% pipeline.
%
%   cellName        (String) name for the context-specific model. It should 
%                   be consistent with the name of the subfolder in which  
%                   the original model is stored. 
%
%   ecModel         Enzyme-constrained model structure without proteomics
%                   constraints.
%   ecModel_batch   Enzyme-constrained model structure with a constrained
%                   total protein pool.
%
% Usage: [ecModel,ecModel_batch] = enhanceGEM_cellLine(cellName)
%
% Ivan Domenzain.      Last edited: 2019-06-10

current = pwd;
cd change_model
model_modified = preprocessModel(model);
% Retrieve kcats & MWs for each rxn in the model from Uniprot database:
cd ../get_enzyme_data
model_data = getEnzymeCodes(model_modified);
%Tries to Match kinetic coefficients to every reaction with a non empty
%grRule
kcats = matchKcats(model_data,'homo sapiens');
cd (current)
mkdir (['../models/' cellName])
cd (current)
%Get ecModel matlab structure
model_data = removeFields(model_data);
cd change_model
ecModel = readKcatData(model_data,kcats);
cd ../limit_proteins
% Constrain total protein pool
Ptotal       = 0.609; %HepG2 total protein content [g prot/gDw]
protCoverage = 0.5;
sigma        = 0.5;
[ecModel_batch,~,~] = constrainEnzymes(ecModel,Ptotal,sigma,protCoverage);
%Save output models:
cd ../../models
ecModel = saveECmodel(ecModel,'COBRA',cellName,version);
ecModel_batch = saveECmodel(ecModel_batch,'COBRA',[cellName '_batch'],version);
cd ../geckomat
end
%--------------------------------------------------------------------------
function model_data = removeFields(model_data)
% Remove unnecessary fields for ecModels that cause conflicts with
% COBRA/RAVEN addRxn functions
%
% Ivan Domenzain.   2018-10-07
%
model = model_data.model;
if isfield(model,'rxnFrom')
    model = rmfield(model,'rxnFrom');
end
if isfield(model,'metFrom')
    model = rmfield(model,'metFrom');
end
if isfield(model,'geneFrom')
    model = rmfield(model,'geneFrom');
end
if isfield(model,'rxnReferences')
    model = rmfield(model,'rxnReferences');
end
if isfield(model,'rxnConfidenceScores')
    model = rmfield(model,'rxnConfidenceScores');
end
if isfield(model,'rxnRecon3DID')
    model = rmfield(model,'rxnRecon3DID');
end
if isfield(model,'rxnRecon3DID')
    model = rmfield(model,'rxnRecon3DID');
end
model_data.model = model;
end