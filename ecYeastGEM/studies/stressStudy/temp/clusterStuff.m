%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clusterStuff(data,name,y_lab,M,N,i,pathways)

%Positions of each element:
if isempty(pathways)
    pos_d  = i;
    pos_hm = (N+i):N:((M-1)*N+i);
else
    pos_d  = 2:N-1;
    pos_hm = [];
    for i = 2:M
        pos_hm = [pos_hm (i*N-2):(i*N-1)];
    end
end

%Remove any row with NaN / with -Inf / with only zeros:
pos  = ((sum(isnan(data),2) > 0) + (sum(isinf(data),2) > 0) + (sum(data,2) == 0));
data = data(~pos,:);

if isempty(pathways)
    %Order from higher sum to lower:
    [~,order] = sort(sum(data,2),'descend');
    data      = data(order,:);
else
    pathways = pathways(~pos);
    %Calculate linkage tree for pathways and create pathway dendrogram:
    tree = linkage(data);
    subplot(M,N,N+1:N:M*N-N+1)
    [H,~,order] = dendrogram(tree,length(pathways),'Orientation','left');
    set(H,'color','k','LineWidth',2)
    ylim([0.5,length(pathways)+0.5])
    axis off
    
    %Reorder data to fit second dendrogram:
    data     = data(order,:);
    pathways = pathways(order);
end

%Calculate linkage tree:
tree = linkage(data');

%Display clusters for PCA:
disp(['Clusters for ' name ':'])
disp(cluster(tree,'maxclust',4)')

%Dendrogram:
subplot(M,N,pos_d)
[H,~,order] = dendrogram(tree);
set(H,'color','k','LineWidth',2)
xlim([0.5,14.5])
axis off

%Heatmap:
% grad = 1000;
% red  = (0:grad)'/grad;
% blue = (grad:-1:0)'/grad;
% map  = [red zeros(grad+1,1) blue];
% colormap(map)
colormap redbluecmap
subplot(M,N,pos_hm)
data = data(end:-1:1,order);
imagesc(data)
colorbar('off')
conds = {'REF','T33','T36','T38','Osm0.2','Osm0.4','Osm0.6','Osm0.8','Osm1.0','Osm1.2','Osm1.3','EtOH20','EtOH40','EtOH60'};
set(gca,'XTick',1:14,'YTick',1:length(pathways),'color','none','YAxisLocation','right')
set(gca,'XTickLabel',conds(order),'YTickLabel',pathways(end:-1:1))
if isempty(pathways)
    ylabel(y_lab,'FontSize',15,'Color','k');
    if N == 3
        set(gca,'FontSize',5)
    end
end

if N > 3
    set(gca,'FontSize',8)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
