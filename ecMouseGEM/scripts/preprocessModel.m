function [model,name,version] = preprocessModel(model,name,version)
%preprocessModel
%
% Performs some preliminary modifications to the metabolic model & 
% retrieves the model's name & version (either by parsing model.id or by
% asking the user to input it), if they were not already defined.
%
% model     A genome-scale model in RAVEN format
% name      The name of the model (alternatively, an empty string)
% version   The version of the model (alternatively, an empty string)
% 
% model     The processed model
% name      The resulting name of the model (if not specified before)
% version   The resulting version of the model (if not specified before)
%
% Usage: [model,name,version] = preprocessModel(model,name,version)

fprintf('Getting genome-scale model ready...')

if nargin< 3
    version = [];
    if nargin <2
        name = [];
    end
end
%correct b vector
model.b = model.b(:,1);

model = removeFields(model);
%find conflicting grRules in model
[~,~,conflicts] = standardizeGrRules(model);
if ~isempty(conflicts)
    warning(sprintf('\nConflicting grRules where found for several reactions (check the above lines). GECKO will ignore these grRules to avoid introduction of potentially wrong enzyme constrains for such cases. It is recommended to check your original GEM and fix these grRules manually with the guidance provided by the function "standardizeGrRules.m" in the RAVEN toolbox.\n'))
end
%delete conflicting grRules to avoid weird enzyme constraints
model.grRules(conflicts) = {''};
%Remove gene rules from pseudoreactions (if any):
pseudoRxns = endsWith(model.rxnNames,' pseudoreaction');
model.grRules(pseudoRxns) = {''};
%standardize grRules
[grRules,rxnGeneMat] = standardizeGrRules(model);
model.grRules        = grRules;
model.rxnGeneMat     = rxnGeneMat;
% Standardizes the metabolites names
model = modifyMetNames(model);

%Delete blocked reactions
%to_remove=model.lb == 0 & model.ub == 0;
%model = removeReactions(model,to_remove,true,true,true);

%Swap direction of only reactions that are defined to only carry negative flux
to_swap=model.lb < 0 & model.ub == 0;
model.S(:,to_swap)=-model.S(:,to_swap);
model.ub(to_swap)=-model.lb(to_swap);
model.lb(to_swap)=0;

%Open all exchange rxns
[~, exchange]      = getExchangeRxns(model);
model.ub(exchange) = +1000;
model.lb(exchange) = -1000;

%Correct rev vector: true if LB < 0 & UB > 0, or if is an exchange reaction:
model.rev = false(size(model.rxns));
for i = 1:length(model.rxns)
    if (model.lb(i) < 0 && model.ub(i) > 0) || sum(model.S(:,i) ~= 0) == 1
        model.rev(i) = true;
    end
end

model = modifyMetNames(model);

if isfield(model,'name')
    name = model.name;
end
if isfield(model,'version')
    version = '0';
end
if isempty(name) && isempty(version) && isfield(model,'id')
    try
        id = strsplit(model.id,'_v');
        if length(id) == 2
            name    = id{1};
            name    = ['ec' upper(name(1)) name(2:end)];
            version = id{2};
        end
    catch
        fprintf('\nNot possible to parse name & version. Input manually\n')
    end
end
while isempty(name)
    name = input('Please enter the desired ecModel name: ','s');
end
while isempty(version)
    version = input('Please enter the model version: ','s');
end
fprintf(' Done!\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = modifyMetNames(model)
% Receives the HMR model .mat structure and modify the metNames field in 
% order to make it compatible with the substrate names of the max_Kcats
% file created with the kinetic data from BRENDA.
%
% Ivan Domenzain. Last edited: 2017-10-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = modifyMetNames(model)
model.OriginalmetNames = model.metNames;
metList                = model.metNames;
for i=1:length(metList)
    met = model.metNames{i};
    % Removes the '_termination code' of each metabolite
    pos = strfind(met,'_');
    if ~isempty(pos)
        met = met(1:pos(end)-1);
    end
    
    if any(strfind(met,'[protein]-'))
        met = replace(met,'[protein]-','');
        
    elseif(any(strfind(met,'-[protein]')))
        met = replace(met,'-[protein]','');
        
    elseif(any(strfind(met,'[protein C terminal]-')))
        met = replace(met,'[protein C terminal]-','');
    end
    model.metNames{i} = met;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = removeFields(model)
if isfield(model,'rxnFrom')
    model = rmfield(model,'rxnFrom');
end
if isfield(model,'inchis')
    model = rmfield(model,'inchis');
end
if isfield(model,'metFormulas')
    model = rmfield(model,'metFormulas');
end
if isfield(model,'metMiriams')
    model = rmfield(model,'metMiriams');
end
if isfield(model,'metCharges')
    model = rmfield(model,'metCharges');
end
if isfield(model,'metFrom')
    model = rmfield(model,'metFrom');
end
if isfield(model,'unconstrained')
    model = rmfield(model,'unconstrained');
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
end






