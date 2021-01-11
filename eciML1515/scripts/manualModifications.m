function [model,modifications] = manualModifications(model)
%manualModifications
%
% Ivan Domenzain.      Last edited: 2019-10-20

%Read manual data:
fID    = fopen('../../databases/manual_data.txt');
data   = textscan(fID,'%s %s %s %s %f %f','delimiter','\t');
rxnIDs = data{1};
kcats  = data{5}.*3600;
SpActs = data{6};
fclose(fID);
modifications{1} = [];
modifications{2} = [];
disp('Improving model with curated data')
for i = 1:length(model.rxns)
    reaction = model.rxnNames{i};
    %Find set of proteins present in rxn:
    S        = full(model.S);
    subs_pos = find(S(:,i) < 0);
    prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
    int_pos  = intersect(subs_pos,prot_pos);
    prot_set = cell(size(int_pos));
    MW_set   = 0;
    for j = 1:length(int_pos)
        met_name    = model.mets{int_pos(j)};
        prot_set{j} = met_name(6:end);
        MW_set      = MW_set + model.MWs(strcmp(model.enzymes,prot_set{j}));
    end
    %Update int_pos:
    S        = full(model.S);
    subs_pos = find(S(:,i) < 0);
    %Get the proteins that are part of the i-th rxn
    prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
    int_pos  = intersect(subs_pos,prot_pos)';
    dataIndx = find(strcmpi(rxnIDs,reaction));
    %For each enzymatic subunit matched to reaction
    for j = 1:length(int_pos)
        enzName = model.mets(int_pos(j));
        if ~isempty(dataIndx)
            kcat  = kcats{dataIndx};
            if ~isempty(kcat)
                newValue = -(kcat)^-1;
                disp(kcat)
            else
                SpAct    = SpActs{dataIndx};
                newValue = -(SpAct*MW_set*0.06)^-1; % [1/h]
                disp(SpAct)
            end
        else
            [newValue,modifications] = curation_topUsedEnz(reaction,enzName,MW_set,modifications);
        end
        %Assign kinetic value in stoichoimetric matrix
        if ~isempty(newValue)
            model.S(int_pos(j),i) = newValue;
        end
    end
    disp(['Improving model with curated data: Ready with rxn #' num2str(i)])
end

%Remove repeated reactions (2017-01-16):
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
    if ~isempty(strfind(rxn_id,'arm_'))
        rxn_code  = rxn_id(5:end);
        k         = 0;
        for j = 1:length(model.rxns)
            if ~isempty(strfind(model.rxns{j},[rxn_code 'No']))
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
modifications = mapModifiedRxns(modifications,model);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation for the model growing on minimal glucose media yielded a list 
% of the top used enzymes (mass-wise), those that were taking more than 5% 
% of the total proteome mass are chosen for manual curation
function [newValue,modifications] = curation_topUsedEnz(reaction,enzName,MW,modifications)
newValue = [];
reaction = string(reaction);
%prot_P0A9D4 (EC2.3.1.30) Value reported for E. coli resulted to be smaller
%in orders of magnitude than the rest of Kcat values for prokaryotes, the
%enzyme was using 14% of the proteome, therefore the highest reported value
%is used instead (salmonella enterica)
if (strcmpi('prot_P0A9D4',enzName)) && contains(reaction,'Serine O-acetyltransferase')
    newValue         = -(200*3600)^-1;
    modifications{1} = [modifications{1}; string('P0A9D4')];
    modifications{2} = [modifications{2}; reaction];
end
%N-acetylglutamate synthase (P0A6C5//E.C.2.3.1.1), enzyme used 9% of the
%proteome, the highest reported value is chosen
if (strcmpi('prot_P0A6C5',enzName)) && contains(reaction,'N-acetylglutamate synthase')
    newValue         = -(0.78*3600)^-1;
    modifications{1} = [modifications{1}; string('P0A6C5')];
    modifications{2} = [modifications{2}; reaction];
end
% (P21151//E.C.2.3.1.16) Enzyme used 7% of proteome mass, Kcat was
% substituted by the highest one
if (strcmpi('prot_P0A6C5',enzName)) && contains(reaction,'Acetyl-CoA C-acyltransferase')
    newValue         = -(14.8*3600)^-1;
    modifications{1} = [modifications{1}; string('P0A6C5')];
    modifications{2} = [modifications{2}; reaction];
end
if (strcmpi('prot_P77399',enzName)) && (contains(reaction,'3-hydroxyacyl-CoA dehydratase') || contains(reaction,'3-hydroxyacyl-CoA dehydrogenase'))
    newValue         = -(6155*MW*60)^-1;
    modifications{1} = [modifications{1}; string('P77399')];
    modifications{2} = [modifications{2}; reaction];
end
%(P05793//E.C.1.1.1.86) There are several values reported  for E. coli and
%this specific E.C. number. The highest reported value is used instead 26
%(1/s) for 3-hydroxypyruvate
if (strcmpi('prot_P05793',enzName)) && (contains(reaction,'Ketol-acid reductoisomerase') || contains(reaction,'2-dehydropantoate 2-reductase'))
    newValue         = -(26*3600)^-1;
    modifications{1} = [modifications{1}; string('P05793')];
    modifications{2} = [modifications{2}; reaction];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modified = mapModifiedRxns(modifications,model)
modified = [];
if ~isempty(modifications)
    for i=1:length(modifications{1})
        rxnIndex = find(strcmp(model.rxnNames,modifications{2}(i)),1);
        str      = {horzcat(modifications{1}{i},'_',num2str(rxnIndex))};
        modified = [modified; str];
    end
end
end
