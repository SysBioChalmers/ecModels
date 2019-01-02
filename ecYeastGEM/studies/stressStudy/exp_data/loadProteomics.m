%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [protIDs,protLevels] = loadProteomics(id,filterProts)

%Read data and transform to correct units:
disp('Reading exp. data files...')
%Merged data (molecules/pgDW):
[~,protIDs]    = xlsread('20170426_merged_proteomic_data.csv',1,'B2:B2319');
[~,conds]      = xlsread('20170426_merged_proteomic_data.csv',1,'C1:AT1');
[protLevels,~] = xlsread('20170426_merged_proteomic_data.csv',1,'C2:AT2319');

%Conversion of units:
protLevels = protLevels*1e12;       %molecules/gDW
protLevels = protLevels/6.02e23;    %mol/gDW
protLevels = protLevels*1e3;        %mmol/gDW

disp('Pre-processing...')
%Replace zeros/negative values by NaN (2015-11-09 after Petri's comment):
disp(['NaN values = ' num2str(sum(sum(isnan(protLevels))))])
disp(['Zero values = ' num2str(sum(sum(protLevels == 0)))])
disp(['Negative values = ' num2str(sum(sum(protLevels < 0)))])
protLevels(protLevels <= 0) = NaN;

%Find specific condition:
pos_cond   = contains(conds,id);
protLevels = protLevels(:,pos_cond);

%Remove any entry with less than 2 measurements:
if filterProts
    nan_pos = false(length(protLevels),1);
    for i = 1:length(protLevels)
        if sum(~isnan(protLevels(i,:))) < 2
            nan_pos(i) = true;
        end
    end
    protLevels(nan_pos,:) = [];
    protIDs(nan_pos,:)    = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
