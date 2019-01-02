%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = manualCuration(model)

to_change = {'r_0118No1','r_1010No1','r_1011No1','r_0904No1','r_2182No1','r_2183No1'};

for i = 1:length(model.rxns)
    if sum(strcmp(to_change,model.rxns{i})) > 0
        %Find set of proteins present in rxn and calculate molecular weight of the set:
        S        = full(model.S);
        subs_pos = find(S(:,i) < 0);
        prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
        int_pos  = intersect(subs_pos,prot_pos);
        MW_set   = 0;
        for j = 1:length(int_pos)
            met_name = model.mets{int_pos(j)};
            MW_set   = MW_set + model.MWs(strcmp(model.enzymes,met_name(6:end)));
        end
        
        %Manual changes:
        for j = 1:length(int_pos)
            % 1. Acetylornithine aminotransferase/transaminase (P18544/EC2.6.1.11):
            % Only kcats available in BRENDA were of S.enterica. Value corrected with
            % s.a. in E.coli  (2016-08-18)
            % Ref: http://www.sciencedirect.com/science/article/pii/0003986179903679
            if strcmp(model.rxns{i},'r_0118No1')
                model.S(int_pos(j),i) = -(153*60*1e3/1e3*MW_set)^-1;    %153 [umol/min/mg]
            end
            
            % 2. Squalene monooxygenase/epoxidase (P32476/EC1.14.14.17=EC1.14.13.132):
            % The kcats available in BRENDA were far too low and none for S.cerevisiae.
            % Value was replaced with the MEDIAN kcat for all EC1.14.14.-, for both 
            % the NADH and NADPH dependant rxns (see stress_problems.xlsx) (2016-08-22)
            if strcmp(model.rxns{i},'r_1010No1') || strcmp(model.rxns{i},'r_1011No1')
                model.S(int_pos(j),i) = -(0.322*3600)^-1;    %0.322 [1/s]
            end

            % 3. Phosphomevalonate kinase (P24521/EC2.7.4.2):
            % The match by GECKO was done to ATP instead of phosphomevalonate because
            % the name in the model was different ((R)-5-phosphomevalonic acid). The
            % value was changed manually to only kcat available for phosphomevalonate
            % in BRENDA (Sus scrofa) (2016-08-19)
            % Ref: http://pubs.acs.org/doi/abs/10.1021/bi00552a004
            if strcmp(model.rxns{i},'r_0904No1')
                model.S(int_pos(j),i) = -(10.2*3600)^-1;    %10.2 [1/s]
            end
            
            % 4. Acyl-CoA desaturase 1 (P21147/EC1.14.19.1): Low kcat values and non of
            % them for yeast. Value corrected with s.a. from Arabidopsis sp. (2017-04-27).
            % Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC55441/
            if strcmp(model.rxns{i},'r_2182No1') || strcmp(model.rxns{i},'r_2183No1')
                model.S(int_pos(j),i) = -(0.8*60*1e3/1e3*MW_set)^-1;    %0.8 [umol/min/mg]
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
