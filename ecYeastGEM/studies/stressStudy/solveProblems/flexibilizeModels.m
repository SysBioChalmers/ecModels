%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [protResults,list] = flexibilizeModels(protResults)

ecModels_wProt = protResults.ecModels_wProt;
list   = cell(size(ecModels_wProt));
[m,n]  = size(ecModels_wProt);
for i = 1:m
    for j = 1:n
        if ~isempty(ecModels_wProt{i,j})
            disp(['Model ' num2str(i) ' - ' num2str(j) ':'])
            vgluc = +100;
            if i == 3
                veth = +100;
            else
                veth = -1;
            end
            model = changeConditions(ecModels_wProt{i,j},vgluc,veth);
            %Solve model problems (so it can grow at least at 0.1 1/h):
            [protResults.ecModels_wProt{i,j},list{i,j}] = flexibilizeModel(model,0.1);
            save('list.mat','list')
            save('../protResults_solved.mat','protResults')
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
