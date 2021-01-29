function [model,GAM] = scaleBioMass(model,Ptot,GAM,scale_comp)
%scaleBioMass
%
% model = scaleBioMass(model,Ptot,GAM,scale_comp)
% 
% Ivan Domenzain.   Last update: 2019-10-20

if nargin < 3
    GAM = [];
end
%Option for changing composition & GAM (=true, default) or only GAM (=false):
if nargin < 4
    scale_comp = true;
end
%Get biomass pseudoreaction ID and biomass components pseudoreactions names
cd ..
parameters = getModelParameters;
cd limit_proteins

[~,Pbase,~,~,~,~] = sumBioMass(model);
%Compute rescaling fractions:
fP = Ptot/Pbase;
%Change compositions:
if scale_comp
	model = rescalePseudoReaction(model,parameters.bio_comp{1},fP);
end
%Change GAM:
xr_pos = strcmp(model.rxns,parameters.bioRxn);
for i = 1:length(model.mets)
    S_ix  = model.S(i,xr_pos);
    isGAM = sum(strcmp({'ATP','ADP','H2O','H+','phosphate'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        GAMpol = 0;
        if isfield(parameters,'pol_cost')
            cost   = parameters.pol_cost;
            GAMpol = Ptot*cost(1) + Ctot*cost(2) + R*cost(3) + D*cost(4);
        end
        model.S(i,xr_pos) = sign(S_ix)*(GAM + GAMpol);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = rescalePseudoReaction(model,metName,f)
rxnName = [metName ' pseudoreaction'];
rxnPos  = strcmp(model.rxnNames,rxnName);
if sum(rxnPos) == 1
    for i = 1:length(model.mets)
        S_ir   = model.S(i,rxnPos);
        isProd = strcmp(model.metNames{i},metName);
        if S_ir ~= 0 && ~isProd
            model.S(i,rxnPos) = f*S_ir;
        end
    end
end
end