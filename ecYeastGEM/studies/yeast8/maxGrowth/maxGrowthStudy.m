% maxGrowthStudy
%   
%   Benjamin J. Sanchez, 2019-02-04
%

%Load model:
ecModel = load('../../../model/ecYeastGEM_batch.mat');
ecModel = ecModel.ecModel_batch;
version = ecModel.id(strfind(ecModel.id,'_v')+1:end);

%Get metabolic model:
git('clone https://github.com/SysBioChalmers/yeast-GEM.git')
cd ./yeast-GEM
git(['checkout ' version])
model = load('./ModelFiles/mat/yeastGEM.mat');
model = ravenCobraWrapper(model.model);
cd ./..

%Load exp data:
fid  = fopen('./maxGrowth_data.csv');
data = textscan(fid,[repmat('%s ',[1,14]) '%s'],'Delimiter',',');
data = [data{1:end}];
fclose(fid);
C_sources = data(1,2:end);
N_sources = data(2:end,1)';
data      = str2double(data(2:end,2:end));

%Block regular consumption (glucose & ammonia):
posGluc = strcmp(ecModel.rxnNames,'D-glucose exchange (reversible)');
posNH4  = strcmp(ecModel.rxnNames,'ammonium exchange (reversible)');
ecModel = setParam(ecModel,'ub',ecModel.rxns(posGluc),0);
ecModel = setParam(ecModel,'ub',ecModel.rxns(posNH4),0);

posGluc = strcmp(model.rxnNames,'D-glucose exchange');
posNH4  = strcmp(model.rxnNames,'ammonium exchange');
model   = setParam(model,'lb',model.rxns(posGluc),0);
model   = setParam(model,'lb',model.rxns(posNH4),0);

%Simulate conditions:
mu_max     = zeros(size(data));
mu_max_pre = zeros(size(data));
fluxes     = cell(size(data));
for i = 1:length(N_sources)
    for j = 1:length(C_sources)
        %Lactate: also open optical isomer:
        if strcmp(C_sources{j},'(R)-lactate')
            vC = [+Inf,+Inf];
        else
            vC = +Inf;
        end
        
        %Run limited model (constraining protein and allowing infinite uptake):
        ecModel_lim = changeMediaCN(ecModel,C_sources{j},N_sources{i},vC,+Inf);
        sol_lim     = solveLP(ecModel_lim,1);
        mu_max(i,j) = -sol_lim.f;
        fluxes{i,j} = sol_lim.x;
        
        %Simulate a normal model:
        posC = strcmp(ecModel.rxnNames,[C_sources{j} ' exchange (reversible)']);
        posN = strcmp(ecModel.rxnNames,[N_sources{i} ' exchange (reversible)']);
        vC   = sol_lim.x(posC);
        vN   = sol_lim.x(posN);
        if strcmp(C_sources{j},'(R)-lactate')
            posC = strcmp(ecModel.rxnNames,'(S)-lactate exchange (reversible)');
            vC   = [vC,sol_lim.x(posC)];
        end
        model_lim       = changeMediaCN(model,C_sources{j},N_sources{i},vC,vN);
        sol_lim         = solveLP(model_lim,1);
        mu_max_pre(i,j) = -sol_lim.f;
        disp(['Growing on ' N_sources{i} ' / ' C_sources{j}])
    end
end

lc = 50;
uc = 500;
rc = 600;
%Fig 1: raw data
figure('position',[lc,lc,rc,uc])
plotMaxGrowth(mu_max,data)

%Figure 2: rescaled data by growth on glucose/ammonia:
data_res = data./data(1,end).*mu_max(1,end);
figure('position',[lc,lc,rc,uc])
plotMaxGrowth(mu_max,data_res)

%Figure 3: rescaled data by growth on glucose:
data_res = data./data(:,end).*mu_max(:,end);
figure('position',[lc,lc,rc,uc])
plotMaxGrowth(mu_max,data_res)

%Figure 4: rescaled data by growth on ammonia:
data_res = data./data(1,:).*mu_max(1,:);
figure('position',[lc,lc,rc,uc])
plotMaxGrowth(mu_max,data_res)

%Figure 5: Breakdown by C source
figure('position',[lc,lc,rc+200,uc])
for i = 1:length(C_sources)
    subplot(3,5,i)
    plotMaxGrowth(mu_max(:,i),data(:,i),C_sources{i},i,[])
end

%Figure 6: Breakdown by N source
figure('position',[lc,lc,rc+200,uc])
for i = 1:length(N_sources)
    subplot(4,6,i)
    plotMaxGrowth(mu_max(i,:),data(i,:),N_sources{i},[],i)
end

%Figure 7: original yeast-GEM:
figure('position',[lc,lc,rc,uc])
plotMaxGrowth(mu_max_pre,data_res)

%Remove the cloned repo:
rmdir('yeast-GEM', 's')
