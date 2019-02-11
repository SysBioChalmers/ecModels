%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model,name,version] = preprocessModel(model,name,version)
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
% Benjamin J. Sanchez. Last edited: 2018-09-01
% Ivan Domenzain.      Last edited: 2019-02-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,name,version] = preprocessModel(model,name,version)
%Convert biomass reaction to a modular type
model = createPoolsForBiomass(model);
%Remove gene rules from pseudoreactions (if any):
for i = 1:length(model.rxns)
    if endsWith(model.rxnNames{i},' pseudoreaction')
        model.grRules{i}      = '';
        model.rxnGeneMat(i,:) = zeros(1,length(model.genes));
    end
end

%Open all exchange rxns
[~, exchange] = getExchangeRxns(model);
model.ub(exchange) = +1000;
model.lb(exchange) = -1000;

%Swap direction of only reactions that are defined to only carry negative flux
to_swap=model.lb < 0 & model.ub == 0;
model.S(:,to_swap)=-model.S(:,to_swap);
model.ub(to_swap)=-model.lb(to_swap);
model.lb(to_swap)=0;

%Delete blocked rxns (LB = UB = 0):
to_remove = logical((model.lb == 0).*(model.ub == 0));
model     = removeReactions(model,model.rxns(to_remove),true,true,true);

%Correct rev vector: true if LB < 0 & UB > 0, or it is an exchange reaction:
model.rev = false(size(model.rxns));
for i = 1:length(model.rxns)
    if (model.lb(i) < 0 && model.ub(i) > 0) || sum(model.S(:,i) ~= 0) == 1
        model.rev(i) = true;
    end
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
        disp('Not possible to parse name & version. Input manually')
    end
end
while isempty(name)
    name = input('Please enter the desired ecModel name: ','s');
end
while isempty(version)
    version = input('Please enter the model version: ','s');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = createPoolsForBiomass(model)
% createPoolsForBiomass
% 
% Create pseudoreactions for carbohydrates, RNA and DNA pools production to
% be incorporated into the biomass pseudoreaction. In this way te biomass
% rxn structure is analogous with the one in yeast8 and biomass composition
% rescaling can be done in a straightforward way.
%
% Ivan Domenzain.  Last edited: 2019-02-08
%

%Biomass components to lump
poolComponents    = cell(1,3);
Pseudometabolites = [{'carbohydrate'},{'RNA'},{'DNA'}];
poolComponents{1} = [{'chitin'},{'trehalose'},{'(1-3)-beta-D-glucan'},{'mannan'}];
poolComponents{2} = [{'AMP'},{'CMP'},{'GMP'},{'UMP'}];
poolComponents{3} = [{'dAMP'},{'dCMP'},{'dGMP'},{'dTMP'}];

%Replace pseudoreaction names 
index                 = find(strcmpi(model.rxnNames,'Amino acid composition'));
model.rxnNames{index} = 'protein pseudoreaction';
index                 = find(strcmpi(model.rxnNames,'Lipids pool'));
model.rxnNames{index} = 'lipid pseudoreaction';

%Add pools as pseudometabolites
metsToAdd.metNames     = Pseudometabolites;
metsToAdd.mets         = Pseudometabolites;
metsToAdd.compartments = [{'c'};{'c'};{'c'}];
metsToAdd.b            = zeros(length(Pseudometabolites),1);
model                  = addMets(model,metsToAdd,false);

%find biomass reaction components
BmPos  = find(strcmpi(model.rxnNames,'Biomass production'));
BmMets = find(model.S(:,BmPos));

for i=1:length(Pseudometabolites)
    PmetIndex    = find(strcmpi(model.metNames,Pseudometabolites{i}));
    components   = poolComponents{i};
    %Identify pool components in metNames
    [~,IA]       = intersect(model.metNames(BmMets),components);
    %Extract original coefficients from biomass reaction
    coefficients = model.S(BmMets(IA),BmPos)';
    %Remove components from biomass pseudoreaction
    model.S(BmMets(IA),BmPos) =  model.S(BmMets(IA),BmPos)*0;
    %add pool to bioRxn
    model.S(PmetIndex,BmPos) = -1;
    %add pseudoreaction for pool production
    rxnName = [Pseudometabolites{i} ' pseudoreaction'];
    rxnsToAdd.rxnNames     = {rxnName};
    rxnsToAdd.rxns         = {rxnName};
    rxnsToAdd.mets         = [model.mets(BmMets(IA)); Pseudometabolites{i}];
    rxnsToAdd.stoichCoeffs = coefficients;
    rxnsToAdd.stoichCoeffs = [rxnsToAdd.stoichCoeffs, 1];
    rxnsToAdd.lb           = 0;
    rxnsToAdd.ub           = 1000;
    rxnsToAdd.c            = 0;
    model = addRxns(model,rxnsToAdd,1,'',false);
end
%model.S = sparse(model.S);
end