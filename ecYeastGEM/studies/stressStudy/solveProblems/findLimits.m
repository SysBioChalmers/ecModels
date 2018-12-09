%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function findLimits(model,option)

if option == 1
    findLimitingEnzyme(model)
elseif option == 2
    findLimitingRxn(model)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function findLimitingEnzyme(model)
base_sol  = optimizeCbModel(model);
prot_rxns = find(~cellfun(@isempty,strfind(model.rxns,'prot_')));
for i = 1:length(prot_rxns)
    if model.ub(prot_rxns(i)) ~= Inf
        new_model = model;
        new_model.ub(prot_rxns(i)) = new_model.ub(prot_rxns(i))*10;
        new_sol   = optimizeCbModel(new_model);
        e         = (new_sol.f - base_sol.f)/base_sol.f;
        if e > 1e-3
            rxn_name = model.rxns{prot_rxns(i)};
            disp(['Protein ' rxn_name(6:11) ': e(mu/P) = ' num2str(e)])
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function findLimitingRxn(model)
base_sol = optimizeCbModel(model);
for i = 1:length(model.rxns)
    new_model = model;
    rxn_mets  = model.mets(full(model.S(:,i)) < 0);
    rxn_prots = find(~cellfun(@isempty,strfind(rxn_mets,'prot_')));
    limiting  = false(size(rxn_prots));
    change    = false;
    for j = 1:length(rxn_prots)
        ex_pos = strcmpi([rxn_mets{rxn_prots(j)} '_exchange'],model.rxns);
        if model.ub(ex_pos) ~= Inf
            new_model.ub(ex_pos) = new_model.ub(ex_pos)*2;
            change               = true;
            limiting(j)          = true;
        end
    end
    if change
        new_sol = optimizeCbModel(new_model);
        e       = (new_sol.f - base_sol.f)/base_sol.f;
        if e > 1e-3
            prots = strrep(rxn_mets(rxn_prots)','prot_','');
            prots = prots(limiting);    %Take out infinite proteins
            disp(['Reaction ' model.rxns{i} ' - Proteins ' strjoin(prots,'-') ...
                  ': e(mu/P) = ' num2str(e)])
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
