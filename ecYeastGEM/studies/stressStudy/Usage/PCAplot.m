%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [all,loadings] = PCAplot(data,name,ecModel)

if nargin == 3
    cd ./../exp_data
    ids  = {'REF','Temp33','Temp36','Temp38','Osmo0.2','Osmo0.4','Osmo0.6', ...
            'Osmo0.8','Osmo1.0','Osmo1.2','Osmo1.3','EtOH20' ,'EtOH40','EtOH60'};
    data = NaN(length(ecModel.enzymes),length(ids));
    for i = 1:length(ids)
        [pIDs,data_i] = loadProteomics(ids{i},false);
        data_i = nanmean(data_i,2);
        for j = 1:length(ecModel.enzymes)
            pos = strcmp(pIDs,ecModel.enzymes{j});
            if sum(pos) == 1
                data(j,i) = data_i(pos);
            end
        end
    end
    cd ./../Usage
    name = 'proteins';
    figure('position', [50,50,500,500]);
end

data = log(data);

%Remove any row with NaN, -Inf or only zeros:
nan_rows  = sum(isnan(data),2) > 0;
inf_rows  = sum(isinf(data),2) > 0;
zero_rows = sum(abs(data),2) == 0;
to_remove = nan_rows + inf_rows + zero_rows > 0;
data      = data(~to_remove,:);

disp(['Number of ' name ' considered in PCA: ' num2str(length(data(:,1)))])

%Return all data (ordered):
[m,n] = size(data);
all   = reshape(data,m*n,1);
all   = sort(all);

%Perform PCA:
[loadings,scores,~,~,explained] = pca(data');

%Reformat loadings matrix:
loadings_ori = zeros(length(to_remove),length(explained));
loadings_ori(~to_remove,:) = loadings;
loadings = loadings_ori;

sizes = [13,13:6:25,13:2:25,13:6:25];
%Plot:
hold on
colors = sampleCVDmap(6);
colors = [[0 0 0];repmat(colors(6,:),3,1);repmat(colors(2,:),7,1);repmat(colors(4,:),3,1)];
for i = length(sizes):-1:1
    plot(scores(i,1),scores(i,2),'o','MarkerEdgeColor','k', ...
         'MarkerFaceColor',colors(i,:),'MarkerSize',sizes(i))
end
sepx  = 3;
sepy  = 1;
minx  = min(min(scores(:,1)));
maxx  = max(max(scores(:,1)));
miny  = min(min(scores(:,2)));
maxy  = max(max(scores(:,2)));
x_min = sepx*floor((minx-(maxx-minx)/5)/sepx);
x_max = sepx*ceil((maxx+(maxx-minx)/10)/sepx);
y_min = sepy*floor((miny-(maxy-miny)/10)/sepy);
y_max = sepy*ceil((maxy+(maxy-miny)/10)/sepy);
setOptions(['PC1: ' num2str(explained(1),3) '% of variation'], ...
            [x_min x_max],x_min:sepx:x_max, ...
            ['PC2: ' num2str(explained(2),3) '% of variation'], ...
            [y_min y_max],y_min:sepy:y_max)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
axis square
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
