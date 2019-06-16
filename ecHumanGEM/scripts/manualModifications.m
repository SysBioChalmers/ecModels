%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = manualModifications(model)
%
% Ivan Domenzain.      Last edited: 2019-06-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = manualModifications(model)

rxnsToremove = {'HMR_9024_REV' 'HMR_9025_REV' 'HMR_9027_REV' 'HMR_9028_REV'...
    'HMR_9029_REV' 'HMR_9030_REV' 'HMR_9031_REV' 'HMR_9035_REV' 'HMR_9036_REV'...
    'HMR_9037_REV' 'HMR_9038_REV' 'HMR_9039_REV' 'HMR_9040_REV' 'HMR_9041_REV'...
    'HMR_9042_REV' 'HMR_9043_REV' 'HMR_9044_REV' 'HMR_9045_REV' 'HMR_9046_REV'...
    'HMR_9048_REV' 'HMR_9049_REV' 'HMR_9051_REV' 'HMR_9052_REV' 'HMR_9053_REV'...
    'HMR_9054_REV' 'HMR_9055_REV' 'humanGrowthOut_REV'};


for i=1:length(rxnsToremove)
    if ~isempty(find(strcmpi(model.rxns,rxnsToremove(i))))
        model = removeReactions(model,rxnsToremove(i),true,true);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%% Other manual changes: %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove repeated reactions (2017-01-16):
rem_rxn = false(size(model.rxns));
for i = 1:length(model.rxns)-1
    for j = i+1:length(model.rxns)
        if isequal(model.S(:,i),model.S(:,j)) && model.lb(i) == model.lb(j) && ...
                model.ub(i) == model.ub(j)
            rem_rxn(j) = true;
            disp(['Removing repeated rxn: ' model.rxns{i} ' & ' model.rxns{j}])
        end
    end
end
model = removeRxns(model,model.rxns(rem_rxn));
% Merge arm reactions to reactions with only one isozyme (2017-01-17):
arm_pos = zeros(size(model.rxns));
p       = 0;
for i = 1:length(model.rxns)
    rxn_id = model.rxns{i};
    if contains(rxn_id,'arm_')
        rxn_code  = rxn_id(5:end);
        k         = 0;
        for j = 1:length(model.rxns)
            if contains(model.rxns{j},[rxn_code 'No'])
                k   = k + 1;
                pos = j;
            end
        end
        if k == 1
            %Condense both reactions in one:
            new_id     = model.rxns{pos};
            new_name   = model.rxnNames{pos};
            stoich     = model.S(:,i) + model.S(:,pos);
            model      = addReaction(model,{new_id,new_name},model.mets,stoich,true,0,1000);
            p          = p + 1;
            arm_pos(p) = i;
            disp(['Merging reactions: ' model.rxns{i} ' & ' model.rxns{pos}])
        end
    end
end
% Remove saved arm reactions:
model = removeRxns(model,model.rxns(arm_pos(1:p)));
% Remove unused enzymes after manual curation (2017-01-16):
rem_enz = false(size(model.enzymes));
for i = 1:length(model.enzymes)
    pos_met = strcmp(model.mets,['prot_' model.enzymes{i}]);
    if sum(model.S(pos_met,:)~=0) == 1
        rem_enz(i) = true;
    end
end
rem_enz = model.enzymes(rem_enz);
for i = 1:length(rem_enz)
    model = deleteProtein(model,rem_enz{i});
    disp(['Removing unused protein: ' rem_enz{i}])
end
end
