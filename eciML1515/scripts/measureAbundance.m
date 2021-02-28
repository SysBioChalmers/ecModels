function [f,count] = measureAbundance(enzymes)
% measureAbundance
%
% Computes the f factor, as a proxy to the mass fraction of proteins accounted 
% for in an ecModel out of the total protein content in cells. An
% integrated quantitative proteomics dataset from the https://pax-db.org/
% database (stored in this toolbox as: 'GECKO/databases/prot_abundance.txt')
% is used as a comparison basis.
% 
% Usage: [f,count] = measureAbundance(enzymes)
%

genes     = {};
abundance = [];
%Read downloaded data of abundance:
%Count number of header lines
fileName = '../../databases/prot_abundance.txt';
if exist(fileName,'file')~= 0
    fID         = fopen('../../databases/prot_abundance.txt');
    data        = textscan(fID,'%s','delimiter','\n');
    if ~isempty(data)
        headerLines = sum(startsWith(data{1},'#'));
        fclose(fID);
        %Read data file, excluding headerlines
        fID         = fopen('../../databases/prot_abundance.txt');
        data        = textscan(fID,'%s %s %f','delimiter','\t','HeaderLines',headerLines);
        genes       = data{2};
        %Remove internal geneIDs modifiers
        genes     = regexprep(genes,'(\d{6}).','');
        abundance = data{3};
        fclose(fID);
    end
else
    fileName = '../../databases/relative_proteomics.txt';
    if exist(fileName,'file')~= 0
        fID       = fopen(fileName);
        data      = textscan(fID,'%s %f','delimiter','\t','HeaderLines',1);
        genes     = data{1};
        abundance = data{2};
    end
end

if ~isempty(genes)
    fprintf('Calculating total abundance...')
    %Load swissprot data:
    data      = load('../../databases/ProtDatabase.mat');
    swissprot = data.swissprot;
    for i = 1:length(swissprot)
        swissprot{i,3} = strsplit(swissprot{i,3},' ');
    end
    
    %Main loop:
    MW_ave  = mean(cell2mat(swissprot(:,5)));
    concs   = zeros(size(genes));
    counter = false(size(genes));
    for i = 1:length(genes)
        MW = MW_ave;
        %Find gene in swissprot database:
        for j = 1:length(swissprot)
            if sum(strcmp(swissprot{j,3},genes{i})) > 0
                MW = swissprot{j,5};	%g/mol
                %Check if uniprot is in model:
                if sum(strcmp(enzymes,swissprot{j,1})) > 0
                    counter(i) = true;
                end
            end
        end
        concs(i) = MW*abundance(i);     %g/mol(tot prot)
        if rem(i,100) == 0 || i == length(genes)
            fprintf('.')
        end
    end
    f     = sum(concs(counter),'omitnan')/sum(concs,'omitnan');
    count = [length(counter);sum(counter,'omitnan')];
else
    disp('prot_abundance file is empty or not available. A default value of f=0.5 is set instead')
    f     = 0.5;
    count = 0;
end
end

