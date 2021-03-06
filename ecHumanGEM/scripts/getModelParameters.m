function parameters = getModelParameters
% getModelParameters
%
%   Set model and organism specific parameters that are used by the
%   ecModel generation pipeline.
%
%   Ivan Domenzain. Last edited: 2020-01-20
%

%Average enzyme saturation factor
parameters.sigma = 0.5;

%Total protein content in the cell [g protein/gDw]
parameters.Ptot = 0.505717472;  %Average across NCI60 cell lines

%Minimum growth rate the model should grow at [1/h]
parameters.gR_exp = 0.020663429; %[g/gDw h]/Average across NCI60 cell lines

%Provide your organism scientific name
parameters.org_name = 'homo sapiens';

%Provide your organism KEGG ID
parameters.keggID = 'hsa';

%The name of the exchange reaction that supplies the model with carbon (rxnNames)
parameters.c_source = 'HMR_9034'; 

%Experimental carbon source uptake (optional)
parameters.c_UptakeExp = 0.641339301; %[mmol/gDw h]/Average across NCI60 cell lines

%Rxn Id for non-growth associated maitenance pseudoreaction
parameters.NGAM = 'biomass_maintenance_Recon3D';

%Compartment name in which the added enzymes should be located
parameters.enzyme_comp = 'Cytosol';
end
