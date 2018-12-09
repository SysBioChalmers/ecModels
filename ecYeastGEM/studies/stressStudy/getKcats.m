%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kcats = getKcats(model)

kcats = zeros(size(model.enzymes));
for i = 1:length(kcats)
    met_pos  = strcmp(model.mets,['prot_' model.enzymes{i}]);
    sub_pos  = model.S(met_pos,:) < 0;
    kcat_set = full(model.S(met_pos,sub_pos));
    kcats(i) = mean(-(kcat_set.^-1)/3600);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
