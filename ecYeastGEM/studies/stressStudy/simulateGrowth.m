%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xS,flux] = simulateGrowth(model,i)

%Find positions of biomass, glucose and protein exchanges:
pos_X    = strcmp(model.rxnNames,'growth');
pos_G    = strcmp(model.rxnNames,'D-glucose exchange (reversible)');
pos_Ein  = strcmp(model.rxnNames,'ethanol exchange (reversible)');
pos_Eout = strcmp(model.rxnNames,'ethanol exchange');
pos_P    = find(strcmp(model.rxnNames,'prot_pool_exchange'));

%Fix biomass, minimize substrate consumption:
model.lb(pos_X) = 0.1;
model.ub(pos_X) = 0.1;
model.lb(pos_G) = 0;
model.ub(pos_G) = 100;
model.c         = zeros(size(model.rxns));
model.c(pos_G)  = -180;                         %6C - MW = 180 g/mol
if i == 3   %Ethanol as substrate
    model.lb(pos_Eout) = 0;
    model.ub(pos_Eout) = 0;
    model.lb(pos_Ein)  = 0;
    model.ub(pos_Ein)  = 100;
    model.c(pos_Ein) = -46;                     %2C - MW = 46 g/mol
end
x = optimizeCbModel(model);

%Afterwards, fix substrates and minimize protein usage:
model.lb(pos_G) = x.x(pos_G);
model.ub(pos_G) = x.x(pos_G)*1.001;
if i == 3   %Ethanol as substrate
    model.lb(pos_Ein) = x.x(pos_Ein);
    model.ub(pos_Ein) = x.x(pos_Ein)*1.001;
end
if isempty(pos_P)   %Minimize fluxes
    model.c = -ones(size(model.rxns));
    x       = optimizeCbModel(model);
else                %Minimize enzymes
    model.c        = zeros(size(model.rxns));
    model.c(pos_P) = -1;
    x              = optimizeCbModel(model);
end

flux = x.x;

%xS:
xS(1) = flux(pos_X);
xS(2) = flux(pos_G);
xS(3) = flux(strcmp(model.rxnNames,'oxygen exchange (reversible)'));
xS(4) = flux(strcmp(model.rxnNames,'carbon dioxide exchange'));
if i == 3
    xS(5) = flux(pos_Ein);
else
    xS(5) = flux(pos_Eout);
end
xS(6) = flux(strcmp(model.rxnNames,'acetate exchange'));
xS(7) = flux(strcmp(model.rxnNames,'glycerol exchange'));
if ~isempty(pos_P)
    xS(8) = flux(pos_P);
else
    xS(8) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
