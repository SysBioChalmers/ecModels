%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,list] = flexibilizeModel(model,required_growth)

%Initialize list of flexibilized enzymes:
list = cell(length(model.enzymes),3);
i    = 0;
K    = 0;

%Calculate current growth:
sol = optimizeCbModel(model);

%Flexibilize until the model grows fast enough:
while sol.f < required_growth
	
	%Find limiting enzymes:
    K = K + 1;
    disp(['  Iteration #' num2str(K) ':'])
    limits = findLimitingEnzyme(model);
	
    if isempty(limits)
        %Find limiting complexes:
        limits = findLimitingRxn(model);
    end
    
    if isempty(limits)
        %Find most used enzyme:
        limits = findMostUsedEnzyme(model);
    end
    
    %Enzyme(s) to add:
    limit = strsplit(limits{1},'-')';
    pos   = (i + 1):(i + length(limit));
    i     = i + length(limit);
    
    %Add new stuff to list:
    list(pos,1) = limit;
    list(pos,2) = getConcentration(model,limit);
    
    %Increase concentration of new stuff:
    model = setConcentration(model,limit,+1000);    %1 mol/gDW
    
    %Re-calculate current growth:
    sol = optimizeCbModel(model);
    disp(['    Current growth: ' num2str(sol.f) ' 1/h'])
end

%Trim unused part:
list = list(1:i,:);

%Minimize usage of flexibilized enzymes with a fixed growth rate:
growth_id = model.rxns(strcmp(model.rxnNames,'growth'));
model_min = changeRxnBounds(model,growth_id,required_growth,'l');
flex_pos  = getEnzymeUsagePos(model_min,list(:,1));
model_min = changeObjective(model_min,model_min.rxns(flex_pos),-1);
sol       = solveLP(model_min,1);

%Re-asign concentrations:
new_concs = max([sol.x(flex_pos) cell2mat(list(:,2))],[],2);
model     = setConcentration(model,list(:,1),new_concs);
list(:,3) = getConcentration(model,list(:,1));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function limits = findLimitingEnzyme(model)

base_sol  = optimizeCbModel(model);
prot_rxns = find(contains(model.rxns,'prot_'));
limits    = cell(0,1);
values    = [];
for i = 1:length(prot_rxns)
    if model.ub(prot_rxns(i)) ~= Inf
        new_model = model;
        new_model.ub(prot_rxns(i)) = new_model.ub(prot_rxns(i))*10;
        new_sol   = optimizeCbModel(new_model);
        e         = (new_sol.f - base_sol.f)/base_sol.f;
        if e > 1e-6
            rxn_name = model.rxns{prot_rxns(i)};
            limits = [limits;rxn_name(6:11)];
            values = [values;e];
            disp(['    Protein ' rxn_name(6:11) ': e(mu/P) = ' num2str(e)])
        end
    end
end
[~,order] = sort(values,'descend');
limits = limits(order);
% limits = limits(~strcmp(limits,'pool_e'));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function limits = findLimitingRxn(model)

base_sol = optimizeCbModel(model);
limits   = cell(0,1);
values   = [];
for i = 1:length(model.rxns)
    new_model = model;
    rxn_mets  = model.mets(full(model.S(:,i)) < 0);
    rxn_prots = find(contains(rxn_mets,'prot_'));
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
        if e > 1e-6
            prots  = strrep(rxn_mets(rxn_prots)','prot_','');
            prots  = prots(limiting);    %Take out infinite proteins
            limits = [limits;strjoin(prots,'-')];
            values = [values;e];
            disp(['    Reaction ' model.rxns{i} ' - Proteins ' strjoin(prots,'-') ...
                  ': e(mu/P) = ' num2str(e)])
        end
    end
end
[~,order] = sort(values,'descend');
limits = limits(order);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function limits = findMostUsedEnzyme(model)

base_sol = solveLP(model,1);
usages   = zeros(length(model.enzymes),2);
for i = 1:length(model.enzymes)
    pos = getEnzymeUsagePos(model,model.enzymes(i));
    if pos > 0
        usages(i,1) = base_sol.x(pos)/model.ub(pos)*100;
        usages(i,2) = base_sol.x(pos)*model.MWs(i)*1000;
    end
end
[~,order] = sortrows(usages,'descend');
limits    = model.enzymes(order);

%Display top 3 most used enzymes:
usages = usages(order,:);
disp(['    Protein ' limits{1} ': Usage = ' num2str(usages(1,1)) '% of ' num2str(usages(1,2)) ' mg/gDW'])
disp(['    Protein ' limits{2} ': Usage = ' num2str(usages(2,1)) '% of ' num2str(usages(2,2)) ' mg/gDW'])
disp(['    Protein ' limits{3} ': Usage = ' num2str(usages(3,1)) '% of ' num2str(usages(3,2)) ' mg/gDW'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function values = getConcentration(model,enzymes)

values = cell(size(enzymes));
for i = 1:length(enzymes)
    pos       = getEnzymeUsagePos(model,enzymes(i));
    values{i} = model.ub(pos);
end
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = setConcentration(model,enzymes,values)

if length(values) == 1
    values = values*ones(size(enzymes));
end

for i = 1:length(enzymes)
    pos   = getEnzymeUsagePos(model,enzymes(i));
    model = changeRxnBounds(model,model.rxns(pos),values(i),'u');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pos = getEnzymeUsagePos(model,enzymes)

pos = zeros(size(enzymes));
for i = 1:length(enzymes)
    try
        pos(i) = find(strcmp(model.rxns,['prot_' enzymes{i} '_exchange']));
    catch
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
