% maxGrowthStudy
%   
%   Benjamin J. Sanchez, 2018-11-22
%

%Load model:
model = load('../../../model/ecYeastGEM_batch.mat');
model = model.ecModel_batch;

%Load exp data:
fid  = fopen('./maxGrowth_data.csv');
data = textscan(fid,[repmat('%s ',[1,14]) '%s'],'Delimiter',',');
data = [data{1:end}];
fclose(fid);
C_sources = data(1,2:end);
N_sources = data(2:end,1)';
data      = str2double(data(2:end,2:end));

%Block regular consumption (glucose & ammonia):
posGluc = strcmp(model.rxnNames,'D-glucose exchange (reversible)');
posNH4  = strcmp(model.rxnNames,'ammonium exchange (reversible)');
posX    = strcmp(model.rxnNames,'growth');
model   = setParam(model,'ub',model.rxns(posGluc),0);
model   = setParam(model,'ub',model.rxns(posNH4),0);

%Simulate conditions:
mu_max = zeros(size(data));
fluxes = cell(size(data));
for i = 1:length(N_sources)
    for j = 1:length(C_sources)
        %Run limited model (constraining protein and allowing infinite uptake):
        model_lim   = changeMediaCN(model,C_sources{j},N_sources{i});
        sol_lim     = solveLP(model_lim,1);
        mu_max(i,j) = sol_lim.x(posX);
        fluxes{i,j} = sol_lim.x;
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
