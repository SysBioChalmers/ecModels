function model = changeMediaCN(model,C_source,N_source)
% changeMedia
%   
%   Usage: model = changeMediaCN(model,C_source,N_source)
%
%   Benjamin J. Sanchez, 2018-11-21
%

%Open substrate uptake rates:
model = openUptake(model,C_source);
model = openUptake(model,N_source);

%Lactate: also open optical isomer:
if strcmp(C_source,'(R)-lactate')
    model = openUptake(model,'(L)-lactate');
end

%Changes to transporters of fructose and mannose: Facilitated instead of
%active (equal to glucose), as concluded from:
%http://onlinelibrary.wiley.com/doi/10.1111/j.1574-6976.1997.tb00346.x/full
if strcmp(C_source,'D-fructose')
    model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1134')) = 0;
    model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1134')) = 0;
elseif strcmp(C_source,'D-mannose')
    model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1139')) = 0;
    model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1139')) = 0;
end

end

function model = openUptake(model,source)

%Open uptake:
pos   = strcmp(model.rxnNames,[source ' exchange (reversible)']);
model = setParam(model,'ub',model.rxns(pos),+Inf);

%Close oposite direction:
pos   = strcmp(model.rxnNames,[source ' exchange']);
model = setParam(model,'ub',model.rxns(pos),0);

end
