%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [path,std] = usageEnrichment(protResults,prot,fluxes,loadings)

delete(get(0,'Children'))   %deletes any present plot

%Apply logarithm to data:
prot.conc = log10(prot.conc);
prot.use  = log10(prot.use);
prot.useP = log10(prot.useP);
prot.flux = log10(fluxes);

[path.names,EPmat,~,RPmat] = getPathways(protResults.ecModels_wMC{1,1});
path.names = strrep(path.names,'and chlorophyll ','');

%Adapt EPmat to fewer data:
EPmat_new = zeros(length(prot.codes),length(path.names));
for i = 1:length(prot.codes)
    pos = strcmp(protResults.ecModels_wMC{1,1}.enzymes,prot.codes{i});
    EPmat_new(i,:) = EPmat(pos,:);
end
EPmat = EPmat_new;

%Find averages of all variables for each pathway:
[path.conc,std.conc] = matchToPathway(prot.conc,path.names,EPmat);
[path.use,std.use]   = matchToPathway(prot.use,path.names,EPmat);
[path.useP,std.useP] = matchToPathway(prot.useP,path.names,EPmat);
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

%Analyze principal components:
path.PC1 = analyzePC(loadings.useP,1,path.names,EPmat);
path.PC2 = analyzePC(loadings.useP,2,path.names,EPmat);

%Over-representation analysis:
path.corrPath  = path.names;
path.corrRatio = zeros(length(path.names),4);
overRep        = false(size(path.names));
for i = 1:length(path.names)
    for j = 1:length(prot.names)
        if sum(strcmp(prot.corrProts,prot.names{j})) > 0
            path.corrRatio(i,1) = path.corrRatio(i,1) + EPmat(j,i);
        else
            path.corrRatio(i,2) = path.corrRatio(i,2) + EPmat(j,i);
        end
    end
    path.corrRatio(i,3) = path.corrRatio(i,1)/(path.corrRatio(i,1) + path.corrRatio(i,2));
    in_corr  = path.corrRatio(i,1);
    in_not   = path.corrRatio(i,2);
    out_corr = length(prot.corrProts) - path.corrRatio(i,1);
    out_not  = length(prot.names) - length(prot.corrProts) - path.corrRatio(i,2);
    path.corrRatio(i,3) = in_corr/(in_corr+in_not);
    if in_corr/(in_corr+in_not) > out_corr/(out_corr+out_not)
        overRep(i) = true;
    end
    Fmatrix = [in_corr,in_not;out_corr,out_not];
    [~,p]   = fishertest(Fmatrix);
    path.corrRatio(i,4) = p;
end
path.corrPath  = strcat(path.corrPath,' (',num2str(sum(path.corrRatio(:,1:2),2)),')');
path.corrPath  = strrep(path.corrPath,'( ','(');
path.corrPath  = path.corrPath(overRep);
path.corrRatio = path.corrRatio(overRep,:);
path.corrPath  = path.corrPath(sum(path.corrRatio(:,1:2),2) > 5);
path.corrRatio = path.corrRatio(sum(path.corrRatio(:,1:2),2) > 5,:);
figure('position',[50,50,1000,250])
hold on
for i = 1:length(path.corrPath)
    score = min(abs(log10(path.corrRatio(i,4)))/log10(20),1);
    barh(i,path.corrRatio(i,3),1,'FaceColor',[score,0,1-score])
    if path.corrRatio(i,4) < 0.01
        pval = ' p < 0.01';
    else
        pval = [' p = ' num2str(round(path.corrRatio(i,4),2))];
    end
    text(path.corrRatio(i,3),i,pval,'HorizontalAlignment','left')
end
limit = length(prot.corrProts)/length(prot.names);
plot([limit limit],[0,length(path.corrPath)]+0.5,'k','LineWidth',2)
setOptions('Fraction of correlated enzymes in pathway',[0,0.4],[],[],[0,length(path.corrPath)]+0.5,1:length(path.corrPath))
set(gca,'YTickLabel',path.corrPath,'FontSize',10)
hold off

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
