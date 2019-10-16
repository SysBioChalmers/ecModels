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
% Ivan Domenzain.      Last edited: 2019-10-15

if nargin< 3
    version = [];
    if nargin <2
        name = [];
    end
end
%Modify some metNames for compatibility with substrate names in the BRENDA
%kinetic data file
model = modifyMetNames(model);
%Changesome reaction names in order to provide compatibility with the
%pipeline
index = find(model.c);
model.rxnNames{index}  = 'biomass pseudoreaction';
index = find(strcmpi(model.rxnNames,'O2 exchange'));
model.rxnNames{index}  = 'oxygen exchange';
index = find(strcmpi(model.rxnNames,'CO2 exchange'));
model.rxnNames{index}  = 'carbon dioxide exchange';
%Introduce biomass pseudometabolite
pseudoMet              = 'biomass';
metsToAdd.metNames     = {pseudoMet};
metsToAdd.mets         = {pseudoMet};
metsToAdd.compartments = {'c'};
metsToAdd.b            = 0;
model                  = addMets(model,metsToAdd,false);
%Incorporate biomass pseudometabolite into biomass pseudoreaction
model.S(end,index)     = 1;
%Introduce exchange reaction for biomass and set as an objective
[model,bioIndex]       = addExchRxn(model,pseudoMet,true);
model.c(bioIndex)      = 1;
%standardize gene-rxn associations
[grRules,rxnGeneMat] = standardizeGrRules(model);
model.grRules        = grRules;
model.rxnGeneMat     = rxnGeneMat;
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
[~, exchange]      = getExchangeRxns(model);
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

%Add name and version as model fields
if isfield(model,'name')
    name = model.name;
end
if isfield(model,'version')
    version = model.version;
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
Pseudometabolites = {'protein'};
components        = [{'Glycine'},{'L-Alanine'},{'L-Leucine'},{'L-Valine'},...
                     {'L-Lysine'},{'L-Arginine'},{'L-Isoleucine'},{'L-Glutamate'},...
                     {'L-Glutamine'},{'L-Threonine'},{'L-Asparagine'},{'L-Aspartate'},...
                     {'L-Proline'},{'L-Serine'},{'L-Phenylalanine'},{'L-Methionine'},...
                     {'L-Tyrosine'},{'L-Histidine'},{'L-Cysteine'},{'L-Tryptophan'}];
%Add pools as pseudometabolites
metsToAdd.metNames     = Pseudometabolites;
metsToAdd.mets         = Pseudometabolites;
metsToAdd.compartments = {'c'};
metsToAdd.b            = zeros(length(Pseudometabolites),1);
model                  = addMets(model,metsToAdd,false);

%find biomass reaction components
BmPos   = find(strcmpi(model.rxnNames,'biomass pseudoreaction'));
BmMets  = find(model.S(:,BmPos));
%Keep the index for "cytosolic biomass only"
Biomass = find(strcmpi(model.metNames,'biomass'),1);
%Add biomass as a product in biomass pseudoreaction
model.S(Biomass,BmPos) = 1;
for i=1:length(Pseudometabolites)
    PmetIndex = find(strcmpi(model.metNames,Pseudometabolites{i}));
    %Identify pool components in metNames
    [~,IA]    = intersect(model.metNames(BmMets),components);
    %Extract original coefficients from biomass reaction
    coefficients = model.S(BmMets(IA),BmPos)';
    %Remove components from biomass pseudoreaction
    model.S(BmMets(IA),BmPos) =  model.S(BmMets(IA),BmPos)*0;
    %add pool to bioRxn
    model.S(PmetIndex,BmPos) = -1;
    %add pseudoreaction for pool production
    rxnName = [Pseudometabolites{i} ' pseudoreaction'];
    rxnsToAdd.rxnNames     = {rxnName};
    rxnsToAdd.rxns         = {[Pseudometabolites{i} '_pool']};
    rxnsToAdd.mets         = [model.mets(BmMets(IA)); Pseudometabolites{i}];
    rxnsToAdd.stoichCoeffs = coefficients;
    rxnsToAdd.stoichCoeffs = [rxnsToAdd.stoichCoeffs, 1];
    rxnsToAdd.lb           = 0;
    rxnsToAdd.ub           = 1000;
    rxnsToAdd.c            = 0;
    model = addRxns(model,rxnsToAdd,1,'',false);
end
model.S = sparse(model.S);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newModel,rxnIndex] = addExchRxn(model,metName,metFlag)

if metFlag
    %Add extracellular metabolite
    metsToAdd.mets          = {[metName '_[e]']};
    metsToAdd.metNames      = {metName};
    metsToAdd.compartments  = {'e'};
    model                   = addMets(model,metsToAdd);
end
%Add transport and exchange reactions
rxnsToAdd.rxns      = {[metName ' transport'];[metName ' exchange']};
rxnsToAdd.rxnNames  = rxnsToAdd.rxns;
rxnsToAdd.equations = {[metName '[c] => ' metName '[e]'];[metName '[e] => ']};
rxnsToAdd.c         = [0 1];
rxnsToAdd.lb        = [0 0];
rxnsToAdd.ub        = [1000 1000];
model.c(:)          = 0;
newModel            = addRxns(model,rxnsToAdd,3);
rxnIndex            = find(strcmpi(model.rxns,[metName ' exchange']));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = modifyMetNames(model)
%Function that modifies metNames for some detected metabolites in the model
%in order to make them compatible with substrate names in the BRENDA
%kinetic data file for the Kcat matching procedure.

index = find(strcmpi(model.metNames,'ADP C10H12N5O10P2'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'ADP';
    end
end

index = find(strcmpi(model.metNames,'ATP C10H12N5O13P3'));
if ~isempty(index)
        for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'ATP';
    end
end

index = find(strcmpi(model.metNames,'CTP C9H12N3O14P3'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'CTP';
    end
end

index = find(strcmpi(model.metNames,'DATP C10H12N5O12P3'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'DATP';
    end
end

index = find(strcmpi(model.metNames,'DCTP C9H12N3O13P3'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'DCTP';
    end
end

index = find(strcmpi(model.metNames,'DGTP C10H12N5O13P3'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'DGTP';
    end
end

index = find(strcmpi(model.metNames,'DTTP C10H13N2O14P3'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'DTTP';
    end
end

index = find(strcmpi(model.metNames,'Fe2+ mitochondria'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'Fe2+';
    end
end

index = find(strcmpi(model.metNames,'GTP C10H12N5O14P3'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'GTP';
    end
end

index = find(strcmpi(model.metNames,'H2O H2O'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'H2O';
    end
end

index = find(strcmpi(model.metNames,'Iron (Fe3+)'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'Fe3+';
    end
end

index = find(strcmpi(model.metNames,'Nicotinamide adenine dinucleotide'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'NAD+';
    end
end

index = find(strcmpi(model.metNames,'Nicotinamide adenine dinucleotide phosphate'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'NADP';
    end
end

index = find(strcmpi(model.metNames,'Nicotinamide adenine dinucleotide phosphate - reduced'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'NADPH';
    end
end

index = find(strcmpi(model.metNames,'Nicotinamide adenine dinucleotide - reduced'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'NADH';
    end
end

index = find(strcmpi(model.metNames,'Flavin adenine dinucleotide reduced'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'FADH2';
    end
end

index = find(strcmpi(model.metNames,'Flavin adenine dinucleotide oxidized'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'FAD';
    end
end

index = find(strcmpi(model.metNames,'Phosphatidylethanolamine (dihexadec-9enoyl, n-C16:1)'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'Phosphatidylethanolamine';
    end
end

index = find(strcmpi(model.metNames,'Phosphatidylethanolamine (dihexadecanoyl, n-C16:0)'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'Phosphatidylethanolamine';
    end
end

index = find(strcmpi(model.metNames,'Protoheme C34H30FeN4O4'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'Protoheme';
    end
end

index = find(strcmpi(model.metNames,'Riboflavin C17H20N4O6'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'Riboflavin';
    end
end

index = find(strcmpi(model.metNames,'Siroheme C42H36FeN4O16'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'Siroheme';
    end
end

index = find(strcmpi(model.metNames,'UTP C9H11N2O15P3'));
if ~isempty(index)
    for i=1:length(index)
        indx = index(i);
        model.metNames{indx} = 'UTP';
    end
end

end