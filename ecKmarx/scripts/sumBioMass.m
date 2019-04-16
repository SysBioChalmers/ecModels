%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,P,C,R,D,L] = sumBioMass(model)
% Calculates breakdown of biomass for the Kmarx model:
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
%        id         MW [g/mol]  class     name
comps = {'s_0955'	89.09       'P'     % A     Alanine         ala
         's_0973'	121.16      'P'     % C     Cysteine        cys
         's_0973'   133.11      'P'     % D     Aspartic acid   asp
         's_0991'   147.13      'P'     % E     Glutamic acid   glu
         's_1032'   165.19      'P'     % F     Phenylalanine   phe
         's_1003'   75.07       'P'     % G     Glycine         gly
         's_1006'   155.15      'P'     % H     Histidine       his
         's_1016'   131.17      'P'     % I     Isoleucine      ile
         's_1025'   146.19      'P'     % K     Lysine          lys
         's_1021'   131.17      'P'     % L     Leucine         leu
         's_1029'   149.21      'P'     % M     Methionine      met
         's_0969'   132.12      'P'     % N     Asparagine      asn
         's_1035'   115.13      'P'     % P     Proline         pro
         's_0999'   146.14      'P'     % Q     Glutamine       gln
         's_0965'   174.2       'P'     % R     Arginine        arg
         's_1039'   105.09      'P'     % S     Serine          ser
         's_1045'   119.12      'P'     % T     Threonine       thr
         's_1056'   117.15      'P'     % V     Valine          val
         's_1048'   204.23      'P'     % W     Tryptophan      trp
         's_1051'   181.19      'P'     % Y     Tyrosine        tyr
         's_0002'	180.16      'C'     % (1->3)-beta-D-glucan
         's_5177'	180.16      'C'     % amylose
         's_0509'   221.21      'C'     % chitin
         's_0773'   180.16      'C'     % glycogen
         's_1107'   180.16      'C'     % mannan
         's_1520'   342.296 	'C'     % trehalose
         's_0423'   347.22      'R'     % AMP
         's_0526'   323.2       'R'     % CMP
         's_0782'   363.22      'R'     % GMP
         's_1545'   324.18      'R'     % UMP
         's_0584'   331.22      'D'     % dAMP
         's_0589'   307.2       'D'     % dCMP
         's_0615'   345.21      'D'     % dGMP
         's_0649'   322.21      'D'     % dTMP
         's_3714'   852.83      'N'     % heme a
         's_1405'   376.36      'N'     % riboflavin
         's_1467'   96.06       'N'     % sulphate
         's_0089'   470.213     'L'     % 1-phosphatidyl-1D-myo-inositol 
         's_0672'   423.65      'L'     % ergosterol ester
         's_0666'   396.65      'L'     % Ergosterol
         's_0694'   259.406     'L'     % fatty acid (average in Yeast; %https://genome.cshlp.org/content/suppl/2003/02/03/13.2.244.DC1/3.pdf)
         's_5125'   226.07      'L'     % phosphatidate
         's_1337'   385.306     'L'     % phosphatidyl-L-serine
         's_1346'   312.085     'L'     % phosphatidylcholine
         's_1351'   269.146     'L'     % phosphatidylethanolamine
         's_1524'   867.22      'L'};   % triglyceride (average in Yeast; %https://genome.cshlp.org/content/suppl/2003/02/03/13.2.244.DC1/3.pdf)
%Get main fractions:
[P,X] = getFraction(model,'protein',0);
[C,X] = getFraction(model,'carbohydrate',X);
[R,X] = getFraction(model,'RNA',X);
[D,X] = getFraction(model,'DNA',X);
[L,X] = getFraction(model,'lipid',X);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,X] = getFraction(model,compType,X)
bioPos     = strcmp(model.rxns,'r_1912');
precursors = model.metNames(find(model.S(:,bioPos)));
precIndex  = find(strcmpi(model.metNames,compType));
F          = abs(model.S(precIndex,bioPos));
X           = X + F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
