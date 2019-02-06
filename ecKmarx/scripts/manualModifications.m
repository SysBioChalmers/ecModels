%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = manualModifications(model)
% 
% Benjamin J. Sanchez. Last edited: 2017-10-29
% Ivan Domenzain.      Last edited: 2018-01-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,modifications] = manualModificationsKma(model)

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
         else
             %%%%%%%%%%%%%%%%%%%% MANUAL CURATION FOR CARBON SOURCES
             [newValue,modifications] = curation_topUsedEnz(reaction,enzName,MW_set,modifications);
             if ~isempty(newValue)
                 model.S(int_pos(j),i) = newValue;
             %else
                 %%%%%%%%%%%%%%% MANUAL CURATION FOR TOP USED ENZYMES:
%                  [newValue,modifications] = curation_carbonSources(reaction,enzName,MW_set,modifications);
%                  if ~isempty(newValue)
%                      model.S(int_pos(j),i) = newValue;
%                  end
             end
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
        
        % 3-hydroxy-3-methylglutaryl coenzyme A reductase (W0TEE1/EC1.1.1.34): 
        % Only kcat available in BRENDA was for Rattus Norvegicus. Value 
        % corrected with max. s.a. in Rattus norvegicus [0.03 umol/min/mg, Mw=226 kDa]
        % from BRENDA (2018-01-27)
          if (strcmpi('prot_W0TEE1',enzName)) %&& (contains(reaction,'(R)-Mevalonate:NADP+ oxidoreductase')))
            % Ratus norvegicus 
            %  newValue         = -(0.03*226000*0.06)^-1;
            % Homarus americanus
             %newValue         = -(60*540000*0.06)^-1; 
            % Streptomyces sp.
             newValue         = -(1.6*100000*0.06)^-1;
             %newValue         = -(60*72000*0.06)^-1;
             modifications{1} = [modifications{1}; string('W0TEE1')];
             modifications{2} = [modifications{2}; reaction];
          end
          
%         % ATP:adenosine 5'-phosphotransferase (W0TF72/EC2.7.1.20): 
%         % No Kcat available for K. marxianus, a value for the same
%         % substrate and for S. cerevisiae was found in BRENDA -> 25.5 1/s
%         % [2018-03-23]
          if (strcmpi('prot_W0TF72',enzName)) && (contains(reaction,'ATP:adenosine 5'))
             newValue         = -(25.5*3600)^-1;
             modifications{1} = [modifications{1}; string('W0TF72')];
             modifications{2} = [modifications{2}; reaction];
          end
%           % dATP:uridine 5'-phosphotransferase W0TCY7//2.7.1.48
%           % dATP:cytidine 5'-phosphotransferase
%           if (strcmpi('prot_W0TCY7',enzName)) %&& (contains(reaction,'ATP:adenosine 5'))
%              newValue         = -(3.2*3600)^-1;
%              modifications{1} = [modifications{1}; string('W0TCY7')];
%              modifications{2} = [modifications{2}; reaction];
%           end
%           %[1.2.1.12] Kcat (Sce & natural substrate) 29 [1/s] is a highly 
%           %growthRate limiting value, S.A. 100 (Sce)
%           %enzIDs = {'prot_P00359','prot_P00360','prot_P00358'};
%           if contains(reaction,'glyceraldehyde-3-phosphate')
%               if (strcmpi('prot_W0T9E0',enzName)) 
%                  newValue = -(234*3600)^-1;
%                  %newValue     = -(29*3600)^-1;
%                  modifications{1} = [modifications{1}; string('W0T9E0')];
%                  modifications{2} = [modifications{2}; reaction];
%               elseif (strcmpi('prot_W0TG93',enzName))
%                  newValue   = -(234*3600)^-1;
%                  newValue     = -(29*3600)^-1;
%                  modifications{1} = [modifications{1}; string('W0TG93')];
%                  modifications{2} = [modifications{2}; reaction];
%               elseif (strcmpi('prot_W0T9I5',enzName))
%                  newValue   = -(234*3600)^-1;
%                  newValue     = -(29*3600)^-1;
%                  modifications{1} = [modifications{1}; string('W0T9I5')];
%                  modifications{2} = [modifications{2}; reaction];
%                
%              end
%          end
%           %[W0T8N1//EC#: 3.3.1.1] S-Adenosyl-L-homocysteine hydrolase (No1)
%           %prev_Kcat:0.79 new_Kcat:6666 gRCC:0.01715 Err:-13.036%
%           % S.A. for Bos taurus used instead
%           if (strcmpi('prot_W0T8N1',enzName))
%             %newValue = (12*237000*60/1000)^-1;
%             modifications{1} = [modifications{1}; string('W0T8N1')];
%             modifications{2} = [modifications{2}; reaction];
%           end
% %           %[W0TDT2//EC#: 2.7.1.40] ATP:pyruvate 2-O-phosphotransferase (No1)
% %           %prev_Kcat:232 new_Kcat:6711.9167 gRCC:0.014911 Err:-9.2003%
% %           %S.A. for Sce used instead
%           if (strcmpi('prot_W0TDT2',enzName))
%             %newValue = (250*3600)^-1;
%             modifications{1} = [modifications{1}; string('W0TDT2')];
%             modifications{2} = [modifications{2}; reaction];
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
    % S-Adenosyl-L-methionine:zymosterol C-methyltransferase [W0TAZ1/EC2.1.1.41]
    % No Kcat reported on BRENDA, the maximum S.A. (E. coli is taken
    % instead)
    if strcmpi('prot_W0TAZ1',enzName) && ...
            (contains(reaction,'zymosterol C-methyltransferase'))
        %newValue      = -(0.01333*3600)^-1;
        %newValue      = -(0.53*172000*60/1000)^-1;
        newValue         = -(0.53*MW_set*60)^-1;
        modifications{1} = [modifications{1}; string('W0TAZ1')];
        modifications{2} = [modifications{2}; reaction];
    end
%     %nicotinate phosphoribosyltransferase [W0TAZ2/EC 6.3.4.21]
%     % S.A. for SCE 2.3 * 45000
%     if strcmpi('prot_W0TAZ2',enzName) %&& (contains(reaction,'zymosterol C-methyltransferase'))
%         newValue         = -(2.3*45000*60/1000)^-1;
%         modifications{1} = [modifications{1}; string('W0TAZ2')];
%         modifications{2} = [modifications{2}; reaction];
%     end
%     % Ketol-acid Reductoisomerase (W0TAV6/EC1.1.1.86)
%     % Substrate name in BRENDA was a synonim as name in model. Changed manually
%     % (2018-03-24).
%     if strcmpi('prot_W0TAV6',enzName) %&& (contains(reaction,'zymosterol C-methyltransferase')))
%         newValue         = -(78.3*3600)^-1;
%         modifications{1} = [modifications{1}; string('W0TAV6')];
%         modifications{2} = [modifications{2}; reaction];
%     end
%     % Amidophosphoribosyltransferase (W0T3E2/EC2.4.2.14)
%     % No reported values for K. marxianus, E. coli was the maximal S.A.
%     % value reported for a microbial organism.
%     % (2018-03-24).
%     if strcmpi('prot_W0T3E2',enzName) && (contains(reaction,'5-phosphoribosylamine:'))
%         newValue         = -(17.2*194000*60/1000)^-1;
%         modifications{1} = [modifications{1}; string('W0T3E2')];
%         modifications{2} = [modifications{2}; reaction];
%     end
%     %Phosphotransferase [W0T3P9/EC 2.7.1.-]
%     % Kcat for Leishmania major used instead
%     if strcmpi('prot_W0T3P9',enzName) 
%         newValue         = -(317*3600)^-1;
%         modifications{1} = [modifications{1}; string('W0T3P9')];
%         modifications{2} = [modifications{2}; reaction];
%     end
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
