%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model,modifications] = manualModificationsGeneral(model)
% 
% Ivan Domenzain.      Last edited: 2019-05-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,modifications] = manualModifications(model)

 modifications{1} = [];
 modifications{2} = [];

%  for i = 1:length(model.rxns)
%      reaction = model.rxnNames{i};
%      %Find set of proteins present in rxn:
%      S        = full(model.S);
%      subs_pos = find(S(:,i) < 0);
%      prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
%      int_pos  = intersect(subs_pos,prot_pos);
%      prot_set = cell(size(int_pos));
%      MW_set   = 0;
%      for j = 1:length(int_pos)
%          met_name    = model.mets{int_pos(j)};
%          prot_set{j} = met_name(6:end);
%          MW_set      = MW_set + model.MWs(strcmp(model.enzymes,prot_set{j}));
%      end
% 
%      %Update int_pos:
%      S        = full(model.S);
%      subs_pos = find(S(:,i) < 0);
%      %Get the proteins that are part of the i-th rxn
%      prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
%      int_pos  = intersect(subs_pos,prot_pos)';
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%  Individual Changes:  %%%%%%%%%%%%%%%%%%%%%%%%
%      for j = 1:length(int_pos)
%          enzName = model.mets(int_pos(j));
%          %%%%%%%%%%%%%%%%%% MANUAL CURATION FOR TOP GROWTH LIMITING ENZYMES:  
%          [newValue,modifications] = curation_growthLimiting(reaction,enzName,MW_set,modifications);
%          if ~isempty(newValue)
%              model.S(int_pos(j),i) = newValue;
%          else
%              %%%%%%%%%%%%%%%%%%%% MANUAL CURATION FOR CARBON SOURCES
%              [newValue,modifications] = curation_topUsedEnz(reaction,enzName,MW_set,modifications);
%              if ~isempty(newValue)
%                  model.S(int_pos(j),i) = newValue;
%              else
%                  %%%%%%%%%%%%%%% MANUAL CURATION FOR TOP USED ENZYMES:
%                  [newValue,modifications] = curation_carbonSources(reaction,enzName,MW_set,modifications);
%                  if ~isempty(newValue)
%                      model.S(int_pos(j),i) = newValue;
%                  end
%              end
%          end          
%       end
%      disp(['Improving model with curated data: Ready with rxn #' num2str(i)])
% end

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
%Change gene rules:
for i = 1:length(model.rules)
    if ~isempty(model.rules{i})
        %Change gene ids:
        model.rules{i} = strrep(model.rules{i},'x(','');
        model.rules{i} = strrep(model.rules{i},')','');
        model.rules{i} = model.genes{str2double(model.rules{i})};
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

modifications = mapModifiedRxns(modifications,model);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the top growth limiting enzymes that were detected by the
% modifyKcats.m script in a preliminary run. 
function [newValue,modifications] = curation_growthLimiting(reaction,enzName,MW_set,modifications)  
        newValue = [];
        reaction = string(reaction);       
        % 3-hydroxy-3-methylglutaryl coenzyme A reductase (W0TEE1/EC1.1.1.34): 
        % Only kcat available in BRENDA was for Rattus Norvegicus. Value 
        % corrected with max. s.a. in Rattus norvegicus [0.03 umol/min/mg, Mw=226 kDa]
        % from BRENDA (2018-01-27)
%           if (strcmpi('prot_W0TEE1',enzName)) && (~isempty(strfind(reaction,'(R)-Mevalonate:NADP+ oxidoreductase (CoA acylating)')))
%              newValue         = -(0.03*226000*0.06)^-1;
%              modifications{1} = [modifications{1}; string('W0TEE1')];
%              modifications{2} = [modifications{2}; reaction];
%           end
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
    if ~isempty(modifications)
        for i=1:length(modifications{1})
            rxnIndex = find(strcmp(model.rxnNames,modifications{2}(i)),1);
            str      = {horzcat(modifications{1}{i},'_',num2str(rxnIndex))};
            modified = [modified; str];
        end
    end
end
