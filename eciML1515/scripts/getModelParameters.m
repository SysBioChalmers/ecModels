function parameters = getModelParameters
% getModelParameters
%
%   Set model and organism specific parameters that are used by the
%   ecModel generation pipeline.
%
%   Ivan Domenzain. Last edited: 2019-07-29
%

%Average enzyme saturation factor
parameters.sigma          = 0.5;
%Total protein content in the cell [g protein/gDw]
parameters.Ptot           = 0.5793;      %Assumed constant
%Minimum growth rate the model should grow at [1/h]
parameters.gR_exp         = 0.58;     %[g/gDw h] 
%Provide your organism scientific name
parameters.org_name       = 'escherichia coli';
%Provide your organism KEGG ID
parameters.keggID         = 'eco';
%The name of the exchange reaction that supplies the model with carbon (rxnNames)
parameters.c_source       = 'D-glucose exchange (reversible)'; 
%Rxn Id for biomass pseudoreaction
parameters.bioRxn         = 'BIOMASS_Ec_iML1515_core_75p37M';
%Compartment name in which the added enzymes should be located
parameters.enzyme_comp    = 'cytosol';
%Rxn names for the most common experimentally measured "exchange" fluxes
parameters.exch_names{1}  = 'biomass exchange';
parameters.exch_names{2}  = 'D-glucose exchange (reversible)';
parameters.exch_names{3}  = 'oxygen exchange (reversible)';
parameters.exch_names{4}  = 'carbon dioxide exchange';
%biomass components pseudoreactions (proteins, carbs and lipids lumped pools)
parameters.bio_comp{1}  = 'protein';
%parameters.bio_comp{2}  = 'carbohydrate';
%parameters.bio_comp{3}  = 'lipid';
%
parameters.NGAM = 1;
end