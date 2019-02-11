%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = manualModifications(model)
%
% Benjamin J. Sanchez. Last edited: 2017-10-29
% Ivan Domenzain.      Last edited: 2018-01-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,modifications] = manualModifications(model)
%Remove reversible growth
model            = removeRxns(model,'y002111_REV');
modifications{1} = [];
modifications{2} = [];

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%  Individual Changes:  %%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:length(int_pos)
        enzName = model.mets(int_pos(j));
        %%%%%%%%%%%%%%%%%% MANUAL CURATION FOR TOP GROWTH LIMITING ENZYMES:
        [newValue,modifications] = curation_growthLimiting(reaction,enzName,MW_set,modifications);
        if ~isempty(newValue)
            model.S(int_pos(j),i) = newValue;
        end
    end
    disp(['Improving model with curated data: Ready with rxn #' num2str(i)])
end
%%%%%%%%%%%%%%%%%%%%%%%%% Other manual changes: %%%%%%%%%%%%%%%%%%%%%%%%%%%
%model = otherChanges(model);
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
%Change gene rules:
if isfield(model,'rules')
    for i = 1:length(model.rules)
        if ~isempty(model.rules{i})
            %Change gene ids:
            model.rules{i} = strrep(model.rules{i},'x(','');
            model.rules{i} = strrep(model.rules{i},')','');
            model.rules{i} = model.genes{str2double(model.rules{i})};
        end
    end
end
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
% Block O2 and glucose production (for avoiding multiple solutions):
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
model.ub(strcmp(model.rxnNames,'D-glucose exchange')) = 0;
% Map the index of the modified Kcat values to the new model (after rxns
% removals).
modifications = mapModifiedRxns(modifications,model);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the top growth limiting enzymes that were detected by the
% modifyKcats.m script in a preliminary run. 
function [newValue,modifications] = curation_growthLimiting(reaction,enzName,MW_set,modifications)  
    newValue = [];
    reaction = string(reaction);
    %*Iteration #1 hexokinase (D-glucose:ATP) [Q6C5S9//EC#: 2.7.1.-]
    if (strcmpi('prot_Q6C5S9',enzName)) && (contains(reaction,'hexokinase (D-glucose:ATP)'))
        % S. cerevisiae S.A. [2018-04-02]
        newValue         = -(120*60000*0.06)^-1; 
        modifications{1} = [modifications{1}; string('Q6C5S9')];
        modifications{2} = [modifications{2}; reaction];       
    end
    %*Iteration #2 lanosterol synthase (No1) [Q6C2X2//EC#:5.4.99.-]
    if (strcmpi('prot_Q6C2X2',enzName)) && (contains(reaction,'lanosterol synthase'))
        % B. taurus S.A. [2018-04-02]
        newValue         = -(1.747*140000*0.06)^-1;   
        modifications{1} = [modifications{1}; string('Q6C2X2')];
        modifications{2} = [modifications{2}; reaction];      
    end
    %*Iteration #3 5'-phosphoribosylformyl glycinamidine synthetase 
    %[Q6BZY1//EC#: 6.3.5.3] prev_Kcat:0.05 new_Kcat:0.35467 gRCC:0.14346 Err:-83.3388%
    if (strcmpi('prot_Q6BZY1',enzName)) && (contains(reaction,'phosphoribosylformyl'))
        % E.coli S.A. [2018-04-02]
        newValue         = -(2.15*141418*0.06)^-1; 
        modifications{1} = [modifications{1}; string('Q6BZY1')];
        modifications{2} = [modifications{2}; reaction];
    end    
    % HMG-CoA reductase (Q6C704/EC1.1.1.34): Only kcat available
    % in BRENDA was for Rattus Norvegicus. Value corrected with max.
    % s.a. in Rattus norvegicus [0.03 umol/min/mg, Mw=226 kDa] from
    % BRENDA (2018-01-27)
    if strcmpi('prot_Q6C704',enzName) && (~isempty(strfind(reaction,'hydroxymethylglutaryl')))
        %(Streptomyces sp. S.A.)
        newValue         = -(1.6*100000*0.06)^-1;
        %newValue         = -(0.03*226000*0.06)^-1;
        modifications{1} = [modifications{1}; string('Q6C704')];
        modifications{2} = [modifications{2}; reaction];
    end
    %[1.2.1.12] 
    % 199 1/s (homo sapiens) used for being the highest non-mutant related
    % value.
    if contains(reaction,'glyceraldehyde-3-phosphate')
        if (strcmpi('prot_Q6CCU7',enzName))
            newValue = -(234*3600)^-1;
            modifications{1} = [modifications{1}; string('Q6CCU7')];
            modifications{2} = [modifications{2}; reaction];
        end
    end
    % 6-phosphofructokinase (P59680/EC2.7.1.11): Kcat for Hsapiens used
    % instead for being the highest value available for the reaction
    % substrates. BRENDA (2018-04-02)
    if strcmpi('prot_P59680',enzName) 
        newValue         = -(357*3600)^-1;
        modifications{1} = [modifications{1}; string('P59680')];
        modifications{2} = [modifications{2}; reaction];
    end
    %Protein:Q6CEM5 // EC#: 2.5.1.6Rxn#:2243 methionine adenosyltransferase 
    %prev_Kcat:8.83e-05 Kcat for R. norvegicus used instead for being the 
    %highest value available for the reaction. 0.583 [1/s]
    %BRENDA (2018-04-18)
    if (strcmpi('prot_Q6CEM5',enzName)) && (contains(reaction,'methionine adenosyltransferase'))
        newValue         = -(0.583*3600)^-1; 
        modifications{1} = [modifications{1}; string('Q6CEM5')];
        modifications{2} = [modifications{2}; reaction];
    end  



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify those kcats involved in extreme misspredictions for growth on 
% several carbon sources. This values were obtained by specific searches on
% the involved pathways for the identification of the ec numbers and then
% its associated Kcat values were gotten from BRENDA.
function [newValue,modifications] = curation_carbonSources(reaction,enzName,MW_set,modifications)
        newValue = [];
        reaction = string(reaction);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After the growth limiting Kcats analysis and the curation for several
% carbon sources, a simulation for the model growing on minimal glucose
% media yielded a list of the top used enzymes (mass-wise), those that were
% taking more than 10% of the total proteome are chosen for manual curation
function [newValue,modifications] = curation_topUsedEnz(reaction,enzName,MW_set,modifications)
    newValue = [];
    reaction = string(reaction);
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = otherChanges(model)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modified = mapModifiedRxns(modifications,model)
    modified = [];
    for i=1:length(modifications{1})
        rxnIndex = find(strcmp(model.rxnNames,modifications{2}(i)),1);
        str      = {horzcat(modifications{1}{i},'_',num2str(rxnIndex))};
        modified = [modified; str];
    end
end