%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,P,C,R,D,L] = sumBioMass(model)
% Calculates breakdown of biomass for the yeast model:
% X -> Biomass fraction without lipids [g/gDW]
% P -> Protein fraction [g/gDW]
% C -> Carbohydrate fraction [g/gDW]
% R -> RNA fraction [g/gDW]
% D -> DNA fraction [g/gDW]
% L -> Lipid fraction [g/gDW]
%
% Ivan Domenzain. Last update: 2019-02-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,P,C,R,D,L] = sumBioMass(model)

%Components of biomass:
%        id             MW [g/mol]  class     name
comps = {'ala__L_c'     89.09       'P'     % A     Alanine         ala
         'cys__L_c'     121.16      'P'     % C     Cysteine        cys
         'CE1787_c'     133.11      'P'     % D     Aspartic acid   asp
         'glu__L_c'     147.13      'P'     % E     Glutamic acid   glu
         'phe__L_c'     165.19      'P'     % F     Phenylalanine   phe
         'gly_c'        75.07       'P'     % G     Glycine         gly
         'his__L_c'     155.15      'P'     % H     Histidine       his
         'ile__L_c'     131.17      'P'     % I     Isoleucine      ile
         'lys__L_c'     146.19      'P'     % K     Lysine          lys
         'leu__L_c'     131.17      'P'     % L     Leucine         leu
         'met__L_c'     149.21      'P'     % M     Methionine      met
         'asn__L_c'     132.12      'P'     % N     Asparagine      asn
         'pro__L_c'     115.13      'P'     % P     Proline         pro
         'gln__L_c'     146.14      'P'     % Q     Glutamine       gln
         'arg__L_c'     174.2       'P'     % R     Arginine        arg
         'ser__L_c'     105.09      'P'     % S     Serine          ser
         'thr__L_c'     119.12      'P'     % T     Threonine       thr
         'val__L_c'     117.15      'P'     % V     Valine          val
         'trp__L_c'     204.23      'P'     % W     Tryptophan      trp
         'tyr__L_c'     181.19      'P'     % Y     Tyrosine        tyr
         '13BDglcn_c'	180.16      'C'     % (1->3)-beta-D-glucan
         'chtn_c'       221.21      'C'     % chitin
         'mannan_c'     180.16      'C'     % mannan
         'tre_c'        342.296 	'C'     % trehalose
         'amp_c'        347.22      'R'     % AMP
         'cmp_c'        323.2       'R'     % CMP
         'gmp_c'        363.22      'R'     % GMP
         'ump_c'        324.18      'R'     % UMP
         'damp_c'       331.22      'D'     % dAMP
         'dcmp_c'       307.2       'D'     % dCMP
         'dgmp_c'       345.21      'D'     % dGMP
         'dtmp_c'       322.21      'D'     % dTMP
         'ergst_c'      396.65      'L'     % Ergosterol
         's_3656_lp'    633.058     'L'     % ergosteryl palmitoleate
         's_3658_lp'    661.112     'L'     % ergosteryl oleate
         'ps_cho_erm'   385.306     'L'     % phosphatidyl-L-serine
         'm1640_lp'     867.22      'L'     % triglyceride (average in Yeast; %https://genome.cshlp.org/content/suppl/2003/02/03/13.2.244.DC1/3.pdf)
         'm1641_lp'     259.406     'L'     % fatty acid (average in Yeast; %https://genome.cshlp.org/content/suppl/2003/02/03/13.2.244.DC1/3.pdf)
         'm1648_erm'    470.213     'L'     % 1-phosphatidyl-1D-myo-inositol
         'm1651_mm'     1466.086    'L'     % cardiolipin
         'm1705_lp'     226.07      'L'     % phosphatidate
         'm1700_lp'     269.146     'L'     % phosphatidylethanolamine
         'm1701_lp'     312.085     'L'};	% phosphatidylcholine
%Get main fractions:
[P,X] = getFraction(model,comps,'P',0);
[C,X] = getFraction(model,comps,'C',X);
[R,X] = getFraction(model,comps,'R',X);
[D,X] = getFraction(model,comps,'D',X);
[L,X] = getFraction(model,comps,'L',X);

%Add up any remaining components:
bioPos = strcmp(model.rxns,'xBIOMASS');
for i = 1:length(model.mets)
    pos = strcmp(comps(:,1),model.mets{i});
    if sum(pos) == 1
        abundance = -model.S(i,bioPos)*comps{pos,2}/1000;
        X         = X + abundance;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,X] = getFraction(model,comps,compType,X)

%Define pseudoreaction name:
rxnName = [compType ' pseudoreaction'];
rxnName = strrep(rxnName,'P','protein');
rxnName = strrep(rxnName,'C','carbohydrate');
rxnName = strrep(rxnName,'R','RNA');
rxnName = strrep(rxnName,'D','DNA');
rxnName = strrep(rxnName,'L','lipid');

%Add up fraction:
fractionPos = strcmp(model.rxnNames,rxnName);
if contains(rxnName,'lipid')
    subs = model.S(:,fractionPos) < 0;        %substrates in pseudo-rxn
    F    = -sum(model.S(subs,fractionPos));   %g/gDW
    F    = F/1000;
else
    comps = comps(strcmp(comps(:,3),compType),:);
    F = 0;
    %Add up all components:
    for i = 1:length(model.mets)
        pos = strcmp(comps(:,1),model.mets{i});
        if sum(pos) == 1
            abundance = -model.S(i,fractionPos)*(comps{pos,2}-18)/1000;
            F         = F + abundance;
        end
    end
end
X = X + F;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
