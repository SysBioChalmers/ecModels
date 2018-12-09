%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sol = testModels(protResults)

%Exp data:
cd ../exp_data
exp_data{1} = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B2:G5');     %Temp
exp_data{2} = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B7:G15');    %Osmotic Stress
exp_data{3} = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B17:G20');   %Ethanol
cd ../solveProblems

models = protResults.ecModels;
[m,n]  = size(models);
sol    = zeros(m,n);
for i = 1:m
    for j = 1:n
        if ~isempty(models{i,j})
            disp(['Condition ' num2str(i) ' - ' num2str(j) ':'])
            %Supply glucose and maybe ethanol:
            vgluc = exp_data{i}(j,1);
            if i == 3
                veth = exp_data{i}(j,4);
            else
                veth = -1;
            end
            model            = manualCuration(models{i,j});
            [model,sol(i,j)] = maxGrowth(model,vgluc,veth);
            if sol(i,j) < 0.1
                findLimits(model,1)
                findLimits(model,2)
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,sol] = maxGrowth(model,vgluc,veth)

%Positions:
pos_X    = strcmp(model.rxnNames,'growth');
pos_G    = strcmp(model.rxnNames,'D-glucose exchange (reversible)');
pos_Ein  = strcmp(model.rxnNames,'ethanol exchange (reversible)');
pos_Eout = strcmp(model.rxnNames,'ethanol exchange');

%Remove any NGAM requirements:
model.lb(strcmp(model.rxnNames,'NGAM')) = 0;

%Fix glucose (and ethanol if it applies) and maximize growth:
model.lb(pos_G) = 0;
model.ub(pos_G) = vgluc;
model.lb(pos_X) = 0;
model.ub(pos_X) = 1000;
model.c         = zeros(size(model.rxns));
model.c(pos_X)  = 1;
if veth >= 0   %Ethanol as substrate
    model.lb(pos_Eout) = 0;
    model.ub(pos_Eout) = 0;
    model.lb(pos_Ein)  = 0;
    model.ub(pos_Ein)  = veth;
end

x   = optimizeCbModel(model);
sol = x.f;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
