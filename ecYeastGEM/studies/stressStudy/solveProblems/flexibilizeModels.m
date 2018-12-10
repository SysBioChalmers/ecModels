%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [protResults,list] = flexibilizeModels(protResults)

%Exp data:
cd ./../exp_data
exp_data{1} = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B2:G5');     %Temp
exp_data{2} = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B7:G15');    %Osmotic Stress
exp_data{3} = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B17:G20');   %Ethanol
cd ./../solveProblems

ecModels_wMC   = protResults.ecModels_wMC;
ecModels_wProt = protResults.ecModels_wProt;
list   = cell(size(ecModels_wProt));
[m,n]  = size(ecModels_wProt);
for i = 1:m
    for j = 1:n
        if ~isempty(ecModels_wProt{i,j})
            disp(['Model ' num2str(i) ' - ' num2str(j) ':'])
            vgluc = exp_data{i}(j,1);
            if i == 3
                veth = exp_data{i}(j,4);
            else
                veth = -1;
            end
            model = changeConditions(ecModels_wMC{i,j},vgluc,veth);
            sol   = optimizeCbModel(model);
            model = changeConditions(ecModels_wProt{i,j},vgluc,veth);
            %Solve model problems (so it can grow at least like wMC):
            [protResults.ecModels_wProt{i,j},list{i,j}] = flexibilizeModel(model,sol.f);
            save('list.mat','list')
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeConditions(model,vgluc,veth)

%Biomass and glucose ids:
X_id    = model.rxns(strcmp(model.rxnNames,'growth'));
G_id    = model.rxns(strcmp(model.rxnNames,'D-glucose exchange (reversible)'));
N_id    = model.rxns(strcmp(model.rxnNames,'non-growth associated maintenance reaction'));
Ein_id  = model.rxns(strcmp(model.rxnNames,'ethanol exchange (reversible)'));
Eout_id = model.rxns(strcmp(model.rxnNames,'ethanol exchange'));

%Fix glucose (and ethanol if applies):
model = changeRxnBounds(model,G_id,0,'l');
model = changeRxnBounds(model,G_id,vgluc,'u');
if veth > 0
    model = changeRxnBounds(model,Ein_id,0,'l');
    model = changeRxnBounds(model,Ein_id,veth,'u');
    model = changeRxnBounds(model,Eout_id,0,'l');
    model = changeRxnBounds(model,Eout_id,0,'u');
end

%Relax NGAM:
model = changeRxnBounds(model,N_id,0,'l');
model = changeRxnBounds(model,N_id,1000,'u');

%Optimize for biomass:
model = changeRxnBounds(model,X_id,0,'l');
model = changeRxnBounds(model,X_id,1000,'u');
model = changeObjective(model,X_id,+1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
