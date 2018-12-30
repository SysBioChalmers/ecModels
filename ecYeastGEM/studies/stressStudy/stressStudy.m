%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function protResults = stressStudy(protResults)

if nargin == 1
    createModels = false;
else
    createModels = true;
end

%Set parameters:
sigma = 0.5;	%Standard value

%Get GECKO:
git clone --depth=1 https://github.com/SysBioChalmers/GECKO.git
cd ./GECKO
folder = fileparts('./*');
addpath(genpath(folder));
GECKOver = git('describe --tags');
cd ./..

%Get metabolic model:
git clone https://github.com/SysBioChalmers/yeast-GEM.git
cd ./yeast-GEM
git fetch --all --tags --prune
git checkout v8.3.1
modelVer = git('describe --tags');
model    = load('./ModelFiles/mat/yeastGEM.mat');
model    = model.model;
cd ./..


%Load and process flux data:
%TODO: switch to csv. [flux_data,samples] = loadFluxData;
cd exp_data
exp_data{1} = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B2:G5');     %Temp
exp_data{2} = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B7:G15');    %Osmotic Stress
exp_data{3} = xlsread('Sce_stress_chemostats_merged.xlsx',1,'B17:G20');   %Ethanol
cd ..

ids = {'REF','Temp33' ,'Temp36' ,'Temp38' ,''       ,''       ,''       ,''       ,''       ;
       'REF','Osmo0.2','Osmo0.4','Osmo0.6','Osmo0.8','Osmo1.0','Osmo1.2','Osmo1.3','Osmo1.4';
       'REF','EtOH20' ,'EtOH40' ,'EtOH60' ,''       ,''       ,''       ,''       ,''      };

%Protein content:
Ptot  = [0.46 0.51 0.52 0.60 NaN  NaN  NaN  NaN  NaN
         0.46 0.63 0.62 0.64 0.62 0.64 0.65 0.64 NaN
         0.46 0.45 0.55 0.52 NaN  NaN  NaN  NaN  NaN];
%TODO: R content + lipids + carbs

%Get enzyme-constrained model:
if createModels
    cd ./GECKO/geckomat
    [ecModel,~] = enhanceGEM(model,'COBRA');
    cd ./../../solveProblems
    ecModel = manualCuration(ecModel);
    cd ./../GECKO/geckomat/limit_proteins
    %model with single enzyme mass constraint:
    Pbase = sumProtein(ecModel);
    [ecModel_general,~,~] = constrainEnzymes(ecModel,Pbase,sigma);
    cd ./../../..
else
    initCobraToolbox
    ecModel_general = protResults.ecModels_wMC{1,1};
end

%Changes to original model:
model       = ravenCobraWrapper(model);
[model,~,~] = preprocessModel(model,'ecYeastGEM','8.0.0');
model       = convertToIrrev(model);
protResults.mModel = model;

%Initialize variables:
M     = length(exp_data);
[m,n] = size(exp_data{2});
protResults.free      = cell(1,M);
protResults.wMC       = cell(1,M);
protResults.wProt     = cell(1,M);
protResults.fluxfree  = cell(1,M);
protResults.fluxwMC   = cell(1,M);
protResults.fluxwProt = cell(1,M);
protResults.protPools = NaN(M,m);
if createModels
    protResults.ecModels_wMC   = cell(M,m);
    protResults.ecModels_wProt = cell(M,m);
end

%For loop for each stress type [i] and each level [j]:
for i = 1:M
    protResults.free{i}      = NaN(m,n+4);
    protResults.wMC{i}       = NaN(m,n+4);
    protResults.wProt{i}     = NaN(m,n+4);
    protResults.fluxfree{i}  = NaN(length(model.rxns),m);
    protResults.fluxwMC{i}   = NaN(length(ecModel_general.rxns),m);
    protResults.fluxwProt{i} = NaN(length(ecModel_general.rxns),m);
    for j = 1:m
        if ~isnan(Ptot(i,j))
            %Changes to original model
            cd ./GECKO/geckomat/limit_proteins
            mModel = scaleBioMass(model,Ptot(i,j));
            
            if createModels
                %Changes to general ecModel:
                ecModel_wMC = scaleBioMass(ecModel_general,Ptot(i,j));
                GAM         = fitGAM(ecModel_wMC);
                cd ./../../..
                
                %Create specific ecModel:
                ecModel_wProt = limitModel(ecModel,ids{i,j},Ptot(i,j),sigma,GAM);
            else
                %Load models:
                ecModel_wMC   = protResults.ecModels_wMC{i,j};
                ecModel_wProt = protResults.ecModels_wProt{i,j};
                cd ./../../..
            end
            
            
            %Fit NGAM to exp data:
            [res_free,flux_free,~]       = iterateNGAM(mModel,exp_data,i,j);
            protResults.free{i}(j,:)     = res_free;
            protResults.fluxfree{i}(:,j) = flux_free;
            
            [res_wMC,flux_wMC,ecModel_wMC] = iterateNGAM(ecModel_wMC,exp_data,i,j);
            protResults.wMC{i}(j,:)        = res_wMC;
            protResults.fluxwMC{i}(:,j)    = flux_wMC;
            
            [res_wProt,flux_wProt,ecModel_wProt] = iterateNGAM(ecModel_wProt,exp_data,i,j);
            protResults.wProt{i}(j,:)            = res_wProt;
            protResults.fluxwProt{i}(:,j)        = flux_wProt;
            
            disp(['Ready with optimization ' num2str(i) ' - ' num2str(j)])
            
            %Store variables:
            protResults.ecModels_wMC{i,j}   = ecModel_wMC;
            protResults.ecModels_wProt{i,j} = ecModel_wProt;
            P_pos = strcmp(ecModel_wProt.rxns,'prot_pool_exchange');
            protResults.protPools(i,j) = ecModel_wProt.ub(P_pos);
            save('protResults.mat','protResults')
        end
    end
end

%Save dependencies:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['GECKO\t' GECKOver '\n']);
fprintf(fid,['yeast-GEM\t' modelVer '\n']);
fclose(fid);

%Remove stuff:
rmpath(genpath(folder));
rmdir('./GECKO','s')
rmdir('./yeast-GEM','s')
delete(get(0,'Children'))   %deletes any present plot

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res,fluxes,model] = iterateNGAM(model,exp_data,i,j)

NGAMspan           = 0:5:50;
[res,~,~]          = optimNGAM(model,exp_data,NGAMspan,i,j);
minNGAM            = max([0,res(end-1)-1]);
NGAMspan_2         = (minNGAM):1:(minNGAM+20);
[res,~,~]          = optimNGAM(model,exp_data,NGAMspan_2,i,j);
minNGAM            = max([0,res(end-1)-1]);
NGAMspan_3         = (minNGAM):0.1:(minNGAM+2);
[res,fluxes,model] = optimNGAM(model,exp_data,NGAMspan_3,i,j);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results,flux,model] = optimNGAM(model,exp_data,NGAM,i,j)

f        = ones(size(NGAM))*1e6;
xS       = cell(size(NGAM));
flux     = cell(size(NGAM));
models   = cell(size(NGAM));
exp_data = exp_data{i}(j,:);

for k = 1:length(NGAM)
    %Change NGAM:
    NGAMpos = strcmp(model.rxnNames,'non-growth associated maintenance reaction');
    model.lb(NGAMpos) = NGAM(k);
    model.ub(NGAMpos) = NGAM(k);
    
    %Simulate growth and calculate f:
    try
        [xS{k},flux{k}] = simulateGrowth(model,i);
        R    = abs(xS{k}(2:end-1)-exp_data);
        f(k) = mean(abs(R));
        disp(['NGAM = ' num2str(NGAM(k)) ' mmol/gDWh -> Error = ' num2str(f(k)) ' mmol/gDWh'])
    catch
        disp('Unfeasible Problem')
        xS{k}   = NaN(1,length(exp_data)+2);
        flux{k} = NaN(length(model.rxns),1);
    end
    models{k} = model;
end

%Choose best:
[~,best] = min(f);
NGAM     = NGAM(best(1));
f        = f(best(1));
xS       = xS{best(1)};
results  = [xS NGAM f];
flux     = flux{best(1)};
model    = models{best(1)};

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
