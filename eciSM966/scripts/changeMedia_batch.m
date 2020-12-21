function [model,pos] = changeMedia_batch(model,c_source,media,flux)
% Kma_changeMedia_batch
%
% function that modifies the ecModel and makes it suitable for batch growth
% simulations on different carbon sources.
%
% model:  An enzyme constrained model
% meadia: Media type ('YEP' for complex, 'MAA' minimal with Aminoacids,
%                          'Min' for minimal media)
% flux:   (Optional) A cell array with measured uptake fluxes in mmol/gDwh
%
% model: a constrained ecModel
%
% usage: [model,pos] = Kma_changeMedia_batch(model,c_source,media,flux)
%
% Ivan Domenzain        2019-02-14

% Give the carbon source (c_source) input variable with the following
% format: c_source  = 'D-glucose exchange (reversible)'

if nargin<3
    media = 'Min';
end
%first block any uptake
[rxnIDs,exchange]  = getExchangeRxns(model);
%Exclude protein pool from exchange reactions list
protIndex = find(contains(model.rxnNames,'prot_'));
exchange  = setdiff(exchange,protIndex);
%First allow any exchange (uptakes and secretions)
model.ub(exchange) = Inf;
%Then block all uptakes
uptakes            = exchange(find(contains(rxnIDs,'_REV')));
model.ub(uptakes)  = 0;
pos = getComponentIndexes(model,c_source);

%Block O2 and glucose production (avoids multiple solutions):
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
model.ub(strcmp(model.rxnNames,'D-glucose exchange')) = 0;
%Find substrate production rxn and block it:
pos_rev = strcmpi(model.rxnNames,c_source(1:strfind(c_source,' (reversible)')-1));
model.ub(pos_rev) = 0;

%The media will define which rxns to fix:
if strcmpi(media,'YEP')
    N = 25;     %Aminoacids + Nucleotides
elseif strcmpi(media,'MAA')
    N = 21;     %Aminoacids
elseif strcmpi(media,'Min')
    N = 1;      %Only the carbon source
end
%UB parameter (manually optimized for glucose on Min+AA):
b = 0.08;
%UB parameter (manually optimized for glucose complex media):
c = 2;
%Define fluxes in case of ec model:
if nargin < 5   %Limited protein    
    if N>1
       flux    = b*ones(1,N);
       if N>21
           flux(22:25) = c;
       end
    end
    flux(1) = Inf;
end
%Fix values as UBs:
for i = 1:N
    model.ub(pos(i)) = flux(i);
end
gIndex = find(model.c);
model.ub(gIndex) = Inf;
%Allow uptake of essential components
model = setParam(model, 'ub', 'r_1727_REV', Inf); % 'ammonium exchange';
model = setParam(model, 'ub', 'r_1723_REV', Inf); % 'water exchange' ;
model = setParam(model, 'ub', 'r_1731_REV', Inf); % 'iron(2+) exchange';
model = setParam(model, 'ub', 'r_1725_REV', Inf); % 'oxygen exchange';
model = setParam(model, 'ub', 'r_1729_REV', 0.1); % 'phosphate exchange';
model = setParam(model, 'ub', 'r_1728_REV', Inf); % 'sulphate exchange';
model = setParam(model, 'ub', 'r_1724_REV', Inf); % 'H+ exchange' ;
model = setParam(model, 'ub', 'r_1772_REV', Inf); % 'Biotin exchange' ;
model = setParam(model, 'ub', 'r_1732_REV', Inf); % 'Nicotinate' ;
model = setParam(model, 'ub', 'r_1730_REV',Inf);  %'myo-inositol'
model = setParam(model, 'ub', 'r_1736_REV',Inf);  %'4-aminobenzoate'
model = setParam(model, 'ub', 'r_1734_REV',Inf);  %'thiamine exchange'
model = setParam(model, 'ub', 'r_1735_REV',Inf);  %'(R)-pantothenate
model = setParam(model, 'ub', 'r_1877_REV',Inf);  %Pyridoxal
model = setParam(model, 'ub', 'r_1771',0);  %Block bicarbonate exchange
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getComponentIndexes(model,c_source)
    pos(1)  = find(strcmpi(model.rxnNames,c_source));
end
