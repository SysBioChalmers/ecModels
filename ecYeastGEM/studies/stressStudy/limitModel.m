%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ecModel_lim,misc_res] = limitModel(ecModel,id,Ptot,sigma,GAM)

%Load and process proteomics data:
cd ./exp_data
[pIDs,data] = loadProteomics(id,true);
cd ./..

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
[ecModel_nostd,~,~] = constrainEnzymes(ecModel,Ptot,sigma,[],GAM,prev_pIDs,data_1);
Pm_nostd = nansum(ecModel_nostd.concs);
Pn_nostd = Ptot - Pm_nostd;

disp('2. MODEL WITH VARIABILITY:')
[ecModel_wstd,~,~] = constrainEnzymes(ecModel,Ptot,sigma,[],GAM,prev_pIDs,data_2);
Pm_wstd = nansum(ecModel_wstd.concs);

disp('3. MODEL WITH VARIABILITY + COMPLEXES FIXED:')
[ecModel_lim,~,~] = constrainEnzymes(ecModel,Ptot,sigma,[],GAM,pIDs,data_3);
Pn_pos = strcmp(ecModel_lim.rxns,'prot_pool_exchange');
Pm_lim = nansum(ecModel_lim.concs);
f_lim  = ecModel_lim.ub(Pn_pos)/Ptot/sigma;
cd ./../../..

disp(['Measured Mass (raw): ' num2str(Pm_nostd) ' g/gDW'])
disp(['Variability introduced: ' num2str(Pm_wstd - Pm_nostd) ' g/gDW'])
disp(['Mass for complexes introduced: ' num2str(Pm_lim - Pm_wstd) ' g/gDW'])
disp(['Total mass in singular constraints: ' num2str(Pm_lim) ' g/gDW'])
disp(['Total mass for general constraint: ' num2str(Pn_nostd) ' g/gDW'])
disp(['fraction f used in general constraint: ' num2str(Ptot*f_lim/Pn_nostd) ' g/g'])
disp(['Corrected mass in general constraint: ' num2str(Ptot*f_lim) ' g/gDW'])
disp(['Total mass in model: ' num2str(Ptot*f_lim + Pm_nostd) ' g/gDW'])
disp(['Mass not in model: ' num2str(Pn_nostd - Ptot*f_lim) ' g/gDW'])

misc_res.Nmatched  = sum(ecModel_nostd.concs > 0);
misc_res.Pmatched  = Pm_nostd;
misc_res.Pstds     = Pm_wstd - Pm_nostd;
misc_res.Pcomplex  = Pm_lim - Pm_wstd;
misc_res.Ppool     = Ptot*f_lim;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
