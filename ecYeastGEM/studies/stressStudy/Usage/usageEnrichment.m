%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [path,std] = usageEnrichment(protResults,prot,fluxes,singleCorr,loadings)

delete(get(0,'Children'))   %deletes any present plot

%Apply logarithm to data:
prot.conc = log10(prot.conc);
prot.useP = log10(prot.useP);
prot.flux = log10(fluxes);

[path.names,EPmat,~,RPmat] = getPathways(protResults.ecModel);

%Adapt EPmat to fewer data:
EPmat_new = zeros(length(prot.codes),length(path.names));
for i = 1:length(prot.codes)
    pos = strcmp(protResults.ecModel.enzymes,prot.codes{i});
    EPmat_new(i,:) = EPmat(pos,:);
end
EPmat = EPmat_new;

%Find averages of all variables for each pathway:
[path.conc,std.conc] = matchToPathway(prot.conc,path.names,EPmat);
[path.useP,std.useP] = matchToPathway(prot.useP,path.names,EPmat);
[path.scor,std.scor] = matchToPathway(singleCorr(:,[2,3]),path.names,EPmat);
[path.flux,std.flux] = matchToPathway(prot.flux,path.names,RPmat);

%Add number of enzymes from pathway:
names      = cell(size(path.names));
path_size  = sum(EPmat,1);
for i = 1:length(path.names)
    names{i}    = [path.names{i} ' (' num2str(path_size(i)) ')'];
end

%Overall top 10 pathways by usage:
useP_means = nanmean(path.useP,2);
useP_stds  = nanstd(path.useP,0,2);     %Std among conditions
filter     = sum(~isnan(path.useP),2) > 3;  %Filter out if half or more of values are NaN.
yellow     = [254 210 24]./255;
topPlot(useP_means(filter),useP_stds(filter),names(filter), 10, ...
      'Pathway usage score',[0 2.5],yellow,true)

%Pathways with usage changing among condition:
corrCond([30,33,36,38],path.useP(:,1:4),names)
corrCond([0:0.2:1.2 1.3],path.useP(:,[1 5:11]),names)
corrCond(0:20:60,path.useP(:,[1 12:14]),names)

%Heriarchical clustering by pathway:
figure('position', [0,0,1000,600])
clusterStuff(path.conc,'conc','Pathways',10,4,1,path.names)
figure('position', [0,0,1000,600])
%Take out outlier:
pos = strcmp(path.names,'');
clusterStuff(path.useP(~pos,:),'useP','Pathways',10,4,1,path.names(~pos))
% clusterStuff(path.useP,'useP','Pathways',10,4,1,path.names)
figure('position', [0,0,1000,600])
clusterStuff(path.flux,'flux','Pathways',10,4,1,path.names)

%Ranking of linear correlations by pathway:
filter = path.scor(:,2) > 0.2;
topPlot(path.scor(filter,1),[],names,10,'Average slope of fit',[],yellow,true)

%Analyze principal components:
path.PC1 = analyzePC(loadings.useP,1,path.names,EPmat);
path.PC2 = analyzePC(loadings.useP,2,path.names,EPmat);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data_path,std_path] = matchToPathway(data,pathways,Pmat)

[m,n]     = size(data);
data_path = NaN(length(pathways),n);
std_path  = NaN(length(pathways),n);
for i = 1:length(pathways)
    for j = 1:n
        pos   = false(1,m);
        for k = 1:m
            if ~isnan(data(k,j)) && ~isinf(data(k,j)) && Pmat(k,i)
                pos(k) = true;
            end
        end
        data_path(i,j) = mean(data(pos,j));
        std_path(i,j)  = std(data(pos,j));
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ranking = analyzePC(loadings,comp,pathways,EPmat)

loadings     = loadings(:,comp);
ranking      = cell(length(pathways),4);
ranking(:,1) = pathways;
for i = 1:length(pathways)
    ranking{i,2} = 0;
    ranking{i,3} = 0;
    for j = 1:length(loadings)
        if EPmat(j,i)
            ranking{i,2} = ranking{i,2} + abs(loadings(j));
            ranking{i,3} = ranking{i,3} + 1;
        end
    end
    if ranking{i,2} > 0
        ranking{i,4} = ranking{i,2}/ranking{i,3};
    else
        ranking{i,4} = 0;
    end
end

%Order results from largest to lowest loading:
[~,order] = sort(cell2mat(ranking(:,2)),'descend');
ranking   = ranking(order,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
