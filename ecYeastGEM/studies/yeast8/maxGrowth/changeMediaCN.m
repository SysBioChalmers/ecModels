function model = changeMediaCN(model,C_source,N_source,C_value,N_value)
% changeMedia
%   
%   Usage: model = changeMediaCN(model,C_source,N_source)
%
%   Benjamin J. Sanchez, 2018-11-21
%

%Open substrate uptake rates:
model = openUptake(model,C_source,C_value);
model = openUptake(model,N_source,N_value);

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

function model = openUptake(model,source,value)

%Lactate: also open optical isomer:
if strcmp(source,'(R)-lactate')
    model = openUptake(model,'(S)-lactate',value(2));
    value = value(1);
end

posB = strcmp(model.rxnNames,[source ' exchange (reversible)']);
posF = strcmp(model.rxnNames,[source ' exchange']);

%Open uptake:
if sum(posB) == 0
    model = setParam(model,'lb',model.rxns(posF),-value);
else
    model = setParam(model,'ub',model.rxns(posB),value);
end

%Close oposite direction:
model = setParam(model,'ub',model.rxns(posF),0);

end
