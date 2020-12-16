function parameters = getModelParameters
% getModelParameters
%
%   Set model and organism specific parameters that are used by the
%   ecModel generation pipeline.
%
%   Ivan Domenzain. Last edited: 2020-09-24
%

%Average enzyme saturation factor
parameters.sigma          = 0.5;
%Total protein content in the cell [g protein/gDw]
parameters.Ptot           = 0.5493;      %Assumed constant
%Minimum growth rate the model should grow at [1/h]
parameters.gR_exp         = 0.278;     %[g/gDw h] 
%Provide your organism scientific name
parameters.org_name       = 'yarrowia lipolytica';
%Provide your organism KEGG ID
parameters.keggID         = 'yli';
%The name of the exchange reaction that supplies the model with carbon (rxnNames)
parameters.c_source       = 'D-glucose exchange (reversible)'; 
%Rxn Id for biomass pseudoreaction
parameters.bioRxn         = 'xBIOMASS';
%Compartment name in which the added enzymes should be located
parameters.enzyme_comp    = 'cytoplasm';
%Rxn names for the most common experimentally measured "exchange" fluxes
parameters.exch_names{1}  = 'growth';
parameters.exch_names{2}  = 'D-glucose exchange (reversible)';
parameters.exch_names{3}  = 'oxygen exchange (reversible)';
parameters.exch_names{4}  = 'carbon dioxide exchange';
%biomass components pseudoreactions (proteins, carbs and lipids lumped pools)
parameters.bio_comp{1}  = 'protein';
parameters.bio_comp{2}  = 'carbohydrate';
parameters.bio_comp{3}  = 'lipid';
%non-growth associated maintenance reaction ID
parameters.NGAM = 'r_4046';
%growth associated maintenance (ATP mmoles/gDw h)
parameters.GAM  = 86.7881;
%Oxphos reaction IDs
parameters.oxPhos{1} = 'y000226';
parameters.oxPhos{2} = 'y000438';
parameters.oxPhos{3} = 'y000439';
parameters.oxPhos{4} = 'y001021';
end