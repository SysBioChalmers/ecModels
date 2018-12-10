%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ecModel_lim = limitModel(ecModel,id,Ptot,sigma,GAM)

%TODO: Load and process proteomics data:
%[protIDs,conds,protLevels] = loadProteomics(Ptot);

%Read data and transform to correct units:
disp('Reading exp. data files...')
%Merged data (molecules/pgDW):
cd exp_data
[~,pIDs]  = xlsread('20170426_merged_proteomic_data.csv',1,'B2:B2319');
[~,conds] = xlsread('20170426_merged_proteomic_data.csv',1,'C1:AT1');
[data,~]  = xlsread('20170426_merged_proteomic_data.csv',1,'C2:AT2319');

%Conversion of units:
data = data*1e12;       %molecules/gDW
data = data/6.02e23;    %mol/gDW
data = data*1e3;        %mmol/gDW
cd ..

disp('Pre-processing...')
%Replace zeros/negative values by NaN (2015-11-09 after Petri's comment):
disp(['NaN values = ' num2str(sum(sum(isnan(data))))])
disp(['Zero values = ' num2str(sum(sum(data == 0)))])
disp(['Negative values = ' num2str(sum(sum(data < 0)))])
data(data <= 0) = NaN;

%Find specific condition and remove any entry with less than 2 measurements:
pos_cond = contains(conds,id);
data     = data(:,pos_cond);
nan_pos  = false(length(data),1);
for i = 1:length(data)
    if sum(~isnan(data(i,:))) < 2
        nan_pos(i) = true;
    end
end
data(nan_pos,:) = [];
pIDs(nan_pos,:) = [];

%Calculate values to use:
data_1    = nanmean(data,2);
dataStds  = nanstd(data,0,2);
data_2    = data_1 + 1.96*dataStds;
data_3    = data_2;
prev_pIDs = pIDs;

%Solve issues in OXPHO complexes (make all measurements proportional to mean):
[data_3,pIDs] = fixComplex('r_1021No1',ecModel,data_3,pIDs);   %Complex II
[data_3,pIDs] = fixComplex('r_0439No1',ecModel,data_3,pIDs);   %Complex III
[data_3,pIDs] = fixComplex('r_0438No1',ecModel,data_3,pIDs);   %Complex IV
[data_3,pIDs] = fixComplex('r_0226No1',ecModel,data_3,pIDs);   %Complex V

%Limit concentrations:
cd ./GECKO/geckomat/limit_proteins
disp('1. MODEL WITH NO EXTRA MASS:')
[eModel_nostd,~,~] = constrainEnzymes(ecModel,Ptot,sigma,[],GAM,prev_pIDs,data_1);
Pm_nostd = nansum(eModel_nostd.concs);
Pn_nostd = Ptot - Pm_nostd;
fm_nostd = Pm_nostd/Ptot;

disp('2. MODEL WITH VARIABILITY:')
[eModel_wstd,~,~] = constrainEnzymes(ecModel,Ptot,sigma,[],GAM,prev_pIDs,data_2);
Pm_wstd = nansum(eModel_wstd.concs);

disp('3. MODEL WITH VARIABILITY + COMPLEXES FIXED:')
[ecModel_lim,~,~] = constrainEnzymes(ecModel,Ptot,sigma,[],GAM,pIDs,data_3);
Pn_pos = strcmp(ecModel_lim.rxns,'prot_pool_exchange');
Pm_lim = nansum(ecModel_lim.concs);
Pn_lim = Ptot - Pm_lim;
f_lim  = ecModel_lim.ub(Pn_pos)/Pn_lim/sigma;
fm_lim = Pm_lim/Ptot;
fn_lim = f_lim*(1 - fm_lim);
cd ./../../..

%Redefine main constraint with the original mass:
f_model = fn_lim/(1-fm_nostd);
ecModel_lim.ub(Pn_pos) = Pn_nostd*f_model*sigma;

disp(['Measured Mass (raw): ' num2str(Pm_nostd) ' g/gDW'])
disp(['Variability introduced: ' num2str(Pm_wstd - Pm_nostd) ' g/gDW'])
disp(['Mass for complexes introduced: ' num2str(Pm_lim - Pm_wstd) ' g/gDW'])
disp(['Total mass in singular constraints: ' num2str(Pm_lim) ' g/gDW'])
disp(['Total mass in general constraint: ' num2str(Pn_nostd) ' g/gDW'])
disp(['fraction f used in general constraint: ' num2str(f_lim) ' g/g'])
disp(['Corrected mass in general constraint: ' num2str(Pn_nostd*f_lim) ' g/gDW'])
disp(['Total mass in model: ' num2str(Pm_lim + Pn_nostd*f_lim) ' g/gDW'])
disp(['Mass not in model: ' num2str(Pn_nostd*(1-f_lim)) ' g/gDW'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
