function [model,pos] = changeMedia_batch(model,c_source,flux)
% changeMedia_batch
%
% function that modifies the ecModel and makes it suitable for batch growth
% simulations on different carbon sources. Script designed for iML1515
% metabolic network of E. coli metabolism.
%
% model     An enzyme constrained model
% c_source (string) Rxn name for the main carbon source uptake reaction
% flux     (Optional) A cell array with measured uptake fluxes in mmol/gDwh
%
% model: a constrained ecModel
%
% usage: [model,pos] = changeMedia_batch(model,c_source,flux)
%
% Ivan Domenzain        2019-10-14

% Provide the carbon source input according to the following format: 
% c_source  = 'D-glucose exchange (reversible)'

if nargin < 3   %Limited protein    
    flux = +1000;
end
pos = find(strcmpi(model.rxnNames,c_source));
%first block any uptake
[rxnIDs,exchange] = getExchangeRxns(model);
%Exclude protein pool from exchange reactions list
protIndex = find(contains(model.rxnNames,'prot_'));
exchange  = setdiff(exchange,protIndex);
%First allow any exchange (uptakes and secretions)
model.ub(exchange) = +1000;
%Then block all uptakes
uptakes            = exchange(find(contains(rxnIDs,'_REV')));
model.ub(uptakes)  = 0;

%Block O2 and glucose production (avoids multiple solutions):
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
%model.ub(strcmp(model.rxnNames,'D-glucose exchange (reversible)')) = 0;
model.ub(strcmp(model.rxnNames,'D-glucose exchange')) = 0;
%Find substrate production rxn and block it:
pos_rev = strcmpi(model.rxnNames,c_source(1:strfind(c_source,' (reversible)')-1));
model.ub(pos_rev) = 0;

%Fix values as UBs:
model.ub(pos)          = flux;
model.ub(find(model.c)) = +1000;
%Allow uptake of essential components
model = setParam(model, 'ub', 'EX_o2_e_REV', +1000); % 'sodium exchange';
%Ions
model = setParam(model, 'ub', 'EX_na1_e_REV', +1000); % 'sodium exchange';
model = setParam(model, 'ub', 'EX_k_e_REV', +1000); % potassium exchange';
model = setParam(model, 'ub', 'EX_zn2_e_REV', +1000); % zinc exchange';
model = setParam(model, 'ub', 'EX_cu_e_REV', +1000); % Cu+ exchange';
model = setParam(model, 'ub', 'EX_cu2_e_REV', +1000); % Cu2+ exchange';
model = setParam(model, 'ub', 'EX_ni2_e_REV', +1000); % Ni2+ exchange';
model = setParam(model, 'ub', 'EX_mn2_e_REV', +1000); % Mn2+ exchange';
model = setParam(model, 'ub', 'EX_mg2_e_REV', +1000); % Mg exchange';
model = setParam(model, 'ub', 'EX_cobalt2_e_REV', +1000); % cobalt exchange';
model = setParam(model, 'ub', 'EX_ca2_e_REV', +1000); % calcium exchange';
model = setParam(model, 'ub', 'EX_mobd_e_REV', +1000); % Molybdate exchange';
model = setParam(model, 'ub', 'EX_fe2_e_REV', +1000); % Fe2+ exchange';
model = setParam(model, 'ub', 'EX_fe3_e_REV', +1000); % Fe3+ exchange';
%Others
model = setParam(model, 'ub', 'EX_pi_e_REV', +1000); % phosphate exchange';
model = setParam(model, 'ub', 'EX_so4_e_REV', +1000); % sulphate exchange';
model = setParam(model, 'ub', 'EX_so3_e_REV', +1000); % sulphite exchange';
model = setParam(model, 'ub', 'EX_so2_e_REV', +1000); % Sulfur dioxide
model = setParam(model, 'ub', 'EX_h2o_e_REV', +1000); % Water exchange
model = setParam(model, 'ub', 'EX_h_e_REV', +1000); % H+ exchange
%Nitrogen source
model = setParam(model, 'ub', 'EX_nh4_e_REV', +1000); % ammonia exchange';
%Inorganic compounds
model = setParam(model, 'ub', 'EX_cl_e_REV', +1000); % chloride exchange';
%Vitamins
model = setParam(model, 'ub', 'EX_thm_e_REV', +1000); % thiamine exchange';
model = setParam(model, 'ub', 'EX_btn_e_REV', +1000); % Biotin exchange';
model = setParam(model, 'ub', 'EX_thym_e_REV', +1000); % Thymine exchange
end