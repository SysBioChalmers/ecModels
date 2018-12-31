
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkExchange(model,sol)

rxns = {};
for i = 1:length(model.rxns)
    pos_met = find(model.S(:,i) ~= 0);
    if sol(i) > 1e-10 && length(pos_met) == 1 && isempty(strfind(model.rxns{i},'prot_'))
        rxn_formula = printRxnFormula(model,model.rxns(i),false,false,true);
        rxns        = [rxns;{sol(i) model.rxns{i} rxn_formula{1}}];
    end
end

[~,order] = sort(cell2mat(rxns(:,1)),'descend');
rxns      = rxns(order,:);
for i = 1:length(order)
    fprintf([sprintf('%6.2e',rxns{i,1}) '\t' rxns{i,2} '\t\t' rxns{i,3} '\n'])
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%