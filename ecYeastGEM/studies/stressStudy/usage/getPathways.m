%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pathways,EPmat,REmat,RPmat] = getPathways(model)

%Construct vector of non repited pathways:
pathways = cell(1000,1);
x        = 1;
for i = 1:length(model.enzymes)
    pathway_i = model.pathways{i};
    pos       = strfind(pathway_i,'sce0');
    for j = 1:length(pos)
        if j == length(pos)
            pos_end = length(pathway_i);
        else
            pos_end = pos(j+1)-2;
        end
        pathways{x} = pathway_i(pos(j):pos_end);
        x           = x + 1;
    end
end
pathways = pathways(1:x-1);
pathways = unique(pathways);
pathways = pathways(1:find(strcmp(pathways,'sce00NNN  None'))-1);

%Create EPmat =true if enzyme i is part of pathway j, =false else:
EPmat = false(length(model.enzymes),length(pathways));
for i = 1:length(model.enzymes)
    for j = 1:length(pathways)
        if ~isempty(strfind(model.pathways{i},pathways{j}))
            EPmat(i,j) = true;
        end
    end
end

%Create REmat =true if rxn i has enzyme j, =false else:
REmat  = false(length(model.rxns),length(model.enzymes));
isProt = ~cellfun(@isempty,strfind(model.mets,'prot_'));
for i = 1:length(model.rxns)
    %Find enzymes corresponding to the rxn:
    enz_pos = find((model.S(:,i) < 0).*isProt);
    for j = 1:length(enz_pos)
        prot_name    = model.mets{enz_pos(j)};
        pos          = strcmp(model.enzymes,prot_name(6:end));
        REmat(i,pos) = true;
    end
end

%Create RPmat =true if rxn i is part of pathway j, =false else:
RPmat  = false(length(model.rxns),length(pathways));
for i = 1:length(model.rxns)
    %Find the enzymes from the reaction:
    enzymes = find(REmat(i,:));
    for j = 1:length(enzymes)
        for k = 1:length(pathways)
            if EPmat(enzymes(j),k)
                RPmat(i,k) = true;
            end
        end
    end
end

%Take out sce code of each pathway:
for i = 1:length(pathways)
    name        = pathways{i};
    pathways{i} = name(11:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
