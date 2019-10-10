% getYeast8FVAfigure
%   
%   Ivan Domenzain, 2019-03-04
%
current = pwd;
c_source = 'r_1714';
%Clone the GECKO repository
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
git ('checkout fix/comparativeFVA')
git pull
cd ..
%Clone the yeast-GEM repository
git('clone https://github.com/SysBioChalmers/yeast-GEM.git')
%Load yeastGEM and constrain it with Y6 media
cd yeast-GEM/ModelFiles/mat
load ('yeastGEM.mat')
%Convert to RAVEN format
model = ravenCobraWrapper(model);
cd ../../ComplementaryScripts/modelCuration
model = minimal_Y6(model);
%Load ecYeastGEM
cd (current)
load ('../../../model/ecYeastGEM_batch.mat')
cd GECKO/geckomat/utilities/FVA
%Run comparative FVA for chemostat conditions
[FVA_chemostat,~,stats_c] = comparativeFVA(model,ecModel_batch,c_source,true,0);
%Run comparative FVA for batch conditions
[FVA_batch,~,stats_b]     = comparativeFVA(model,ecModel_batch,c_source,false,0);
%convert zeros to low values to get a visually nice plot
FVA_chemostat = filterDistributions(FVA_chemostat,1E-10);
FVA_batch     = filterDistributions(FVA_batch,1E-10);
distributions = {FVA_chemostat{1}, FVA_chemostat{2},FVA_batch{1},FVA_batch{2}};
stats         = {stats_c,stats_b}; 
legends       = {'model-chemostat', 'ecModel-chemostat','model-batch', 'ecModel-batch'};
titleStr      = 'Flux variability cumulative distribution';
[~, ~]        = plotCumDist(distributions,legends,titleStr);
%Remove the cloned repos:
cd (current)
rmdir('GECKO', 's')
rmdir('yeast-GEM', 's')

function newDist = filterDistributions(distribution,treshold)
newDist = [];
for i=1:length(distribution)
    dist = distribution{i};
    dist(dist<treshold) =treshold;
    dist(find(dist==1E-10,1)) = treshold/10;
    newDist = [newDist,{dist}];
end
end