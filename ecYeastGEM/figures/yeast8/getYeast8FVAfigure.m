function [distributions, stats] = getYeast8FVAfigure(model,ecModel)
c_source = 'D-glucose exchange';
[FVA_chemostat,~,stats_c] = comparativeFVA(model,ecModel,c_source,true,0,'oxygen');
[FVA_batch,~,stats_b] = comparativeFVA(model,ecModel,c_source,false,0,'oxygen');
distributions = {FVA_chemostat{1}, FVA_chemostat{2},FVA_batch{1},FVA_batch{2}};
stats         = {stats_c,stats_b}; 
legends       = {'model-chemostat', 'ecModel-chemostat','model-batch', 'ecModel-batch'};
titleStr      = 'Flux variability cumulative distribution';
[~, ~]        = plotCumDist(distributions,legends,titleStr);
end