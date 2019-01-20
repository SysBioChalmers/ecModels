%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [protIDs,protLevels] = loadProteomics(id,filterProts,Ptot)

%Read data and transform to correct units:
disp('Reading exp. data files...')
%Merged data:
[~,protIDs]    = xlsread('merged_proteomic_data.csv',1,'B2:B6000');
[~,conds]      = xlsread('merged_proteomic_data.csv',1,'E1:BZ1');
[protLevels,~] = xlsread('merged_proteomic_data.csv',1,'E2:BZ6000');    %g/g detected
[MWs,~]        = xlsread('merged_proteomic_data.csv',1,'D2:D6000');     %kDa = g/mmol

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
    MWs(nan_pos,:)        = [];
    %Rescale again to 1:
    protLevels = protLevels./nansum(protLevels);
end

%Conversion of units:
cd ./../GECKO/geckomat/limit_proteins
[f,~] = measureAbundance(protIDs);  %Fraction of proteins detected
cd ./../../../exp_data
protLevels = protLevels*f;          %g/g protein
protLevels = protLevels*Ptot;       %g/gDW
protLevels = protLevels./MWs;       %mmol/gDW

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
