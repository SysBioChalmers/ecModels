%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = scaleBioMass(model,Ptot,GAM,scale_comp)
% 
% Ivan Domenzain.   Last update: 2019-06-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,GAM] = scaleBioMass(model,Ptot,GAM,scale_comp)
if nargin < 3
    GAM = [];
end
%Option for changing composition & GAM (=true, default) or only GAM (=false):
if nargin < 4
    scale_comp = true;
end
[~,Pbase,~,~,~,~] = sumBioMass(model);
%Compute rescaling fractions:
fP = Ptot/Pbase;
%Change compositions:
if scale_comp
    model = rescalePseudoReaction(model,'protein',fP);
end
%Fit GAM if not available:
if isempty(GAM)
    xr_pos  = strcmp(model.rxnNames,'biomass pseudoreaction');
    S_ix    = find(model.S(:,xr_pos));
    bioMets = model.metNames(S_ix);
    index   = S_ix(strcmpi(bioMets,'ATP'));
    if ~isempty(index)
        GAM = abs(model.S(index,xr_pos));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
