function exchModel = changeMedia_batch(model,cSource,irrev,measuredMets,fluxes)
% changeMedia_batch
%
% Set a Hams culture medium for a humanGEM based model. This function works
% for either standard or enzyme constrained-GEMs
%
%   model           An ihuman-based GEM
%   irrev           (logical) TRUE if model comes in an irreversible format.
%                   Default = false
%   measuredMets    (string) metNames for measured compounds. Optional
%   fluxes          (vector) Measured fluxes [mmol/gDw h]. Optional
%                   NOTE: NEGATIVE = CONSUME
%                         POSITIVE = PRODUCE
%
%   exchModel       (struct) Model with Ham's media constraints
%
% Ivan Domenzain.      Last edited: 2020-11-20

if nargin<5
	fluxes = [];
    if nargin<4
    	measuredMets = [];
        if nargin<3 || isempty(irrev)
        	irrev = true;
        end
    end
end
%Remove unconstrained field, if available
if isfield(model,'unconstrained')
    model = rmfield(model,'unconstrained');
end
%Check if boundary metabolites are present, if so then remove them
boundaryIndx = find(strcmpi(model.compNames,'Boundary'));
if ~isempty(boundaryIndx)
    boundary = find(model.metComps==boundaryIndx);
    model    = removeMets(model,boundary,false,false,false,true);
end
%Ham's media composition
mediaComps ={'glucose'
             'arginine'
             'histidine'
             'lysine'
             'methionine'
             'phenylalanine'
             'tryptophan'
             'tyrosine'
             'alanine'
             'glycine'
             'serine'
             'threonine'
             'aspartate'
             'glutamate'
             'asparagine'
             'glutamine'
             'isoleucine'
             'leucine'
             'proline'
             'valine'
             'cysteine'
             'thiamin'
             'hypoxanthine'
             'folate'
             'biotin'
             'pantothenate'
             'choline'
             'inositol'
             'nicotinamide'
             'pyridoxine'
             'riboflavin'
             'thymidine'
             'aquacob(III)alamin'
             'lipoic acid'
             'sulfate'
             'linoleate'
             'linolenate'
             'O2'
             'H2O'
             'retinoate'
             'Fe2+'
             'Pi'
             'alpha-tocopherol'
             'gamma-tocopherol'};
             
%Default flux bounds [LB, UB]
fluxBounds = [-ones(length(mediaComps),1), ones(length(mediaComps),1)]*1000;
%Check if provided mets are part of media's formulation
if ~isempty(measuredMets)
    %Modify fluxBounds with the correspondent provided flux measurements
    [iA,iB] = ismember(measuredMets,mediaComps);
    fluxBounds(iB(iA),:) = fluxes(iA).*ones(sum(iA),2);  % force flux to be equal to measured value (LB = UB = measured val)
    if any(~iA)
        %If measured mets are not in media formulation, then add them
        mediaComps = [mediaComps; measuredMets(~iA)];
        fluxBounds = [fluxBounds; fluxes(~iA).*ones(sum(~iA),2)];
    end
end

if ~irrev
    modelStr = 'model';
    %Set uptake fluxes for media mets
    [exchModel,unusedMets] = setExchangeBounds(model,mediaComps,fluxBounds(:,1),fluxBounds(:,2),true);
else
    modelStr = 'ecModel';
    [exchRxns,exchIndxs] = getExchangeRxns(model);
    %Exclude protein pool exchange
    exchIndxs = exchIndxs(1:end-1);
    exchRxns  = exchRxns(1:end-1);
    %Differentiate between uptakes and production reactions
    uptkIndxs = exchIndxs(find(contains(exchRxns,'_REV')));
    prodIndxs = exchIndxs(find(~contains(exchRxns,'_REV')));
    %Open all production reactions
    exchModel = setParam(model,'ub',prodIndxs,1000);
    %close all uptakes
    exchModel = setParam(exchModel,'ub',uptkIndxs,0);
    %Open uptake of media components one by one
    unusedMets = [];
    for i=1:length(mediaComps)
        %Get metabolite indx
        metIndx = getIndexes(model,strcat(mediaComps{i},'[e]'),'metcomps');
        %Get rxns for metabolite
        metRxns = find(model.S(metIndx,:));
        %Get the uptake reaction for the metabolite
        metExchRxn = intersect(metRxns,uptkIndxs);
        if isempty(metExchRxn)
            unusedMets = [unusedMets; mediaComps(i)];
        else
            %Open it!
            if fluxBounds(i,1) < 0
                % if the metabolite is being consumed
                exchModel.ub(metExchRxn) = abs(fluxBounds(i,1));
            else
                exchModel.ub(metExchRxn) = 0;
            end
            %If the metabolite was measured, also set its production rate
            metExchRxn = intersect(metRxns,prodIndxs);
            if fluxBounds(i,2) > 0
                exchModel.ub(metExchRxn) = abs(fluxBounds(i,2));
            else
                exchModel.ub(metExchRxn) = 0;
            end
        end
    end
end

% report unused metabolites
if ~isempty(unusedMets)
    fprintf('WARNING: The following metabolites are either not in the model or do not have exchange reactions:\n');
    fprintf('\t%s\n',unusedMets{:});
end

%Check if model is feasible
sol = solveLP(exchModel); 
if ~isempty(sol.x)
	disp(['Constrained ' modelStr ' is feasible'])
else
    disp(['Constrained ' modelStr ' is unfeasible'])
	exchModel = [];
end
end