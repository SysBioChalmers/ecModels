function filesGeneIDs
current = pwd;
%Load ecYeast model
cd ../../model
load('ecYeastGEM_batch.mat')
MediaTypes = {'Min' 'MAA' 'YEP'};
cd (current)
for media = MediaTypes
     fileName   = ['KcatSensitivities_' media{1}, '.txt'];
     cd ([current '/results'])
     dataTable = readtable(fileName);
     [~,n] = size(dataTable);
     [~,indexes]  = ismember(dataTable.Protein,ecModel_batch.enzymes);
     genes = ecModel_batch.enzGenes(indexes);
     newTable = {};
     for i=2:n
        newTable = [newTable,dataTable(:,i)];
     end
     newTable = [genes, newTable];
     fileName   = ['KcatSensitivities_' media{1} '_genes', '.txt'];
     writetable(newTable, fileName,'Delimiter','\t','QuoteStrings',false,'WriteRowNames',false);
end
end