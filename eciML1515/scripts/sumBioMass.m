function [X,P,C,L,D,R] = sumBioMass(model)
%sumBioMass
%
% Calculates breakdown of biomass for the iML1515 model:
% X -> Biomass fraction without lipids [g/gDW]
% P -> Protein fraction [g/gDW]
%
% Usage: [X,P] = sumBioMass(model)
%
% Ivan Domenzain. Last update: 2019-10-20


%Components of biomass:
%        id         MW [g/mol]  class     name
comps = {'ala__L_c'	  89.09      'P'     % A     Alanine         ala
         'cys__L_c'	 121.16      'P'     % C     Cysteine        cys
         'asp__L_c'  133.11      'P'     % D     Aspartic acid   asp
         'glu__L_c'  147.13      'P'     % E     Glutamic acid   glu
         'phe__L_c'  165.19      'P'     % F     Phenylalanine   phe
         'gly_c'      75.07      'P'     % G     Glycine         gly
         'his__L_c'  155.15      'P'     % H     Histidine       his
         'ile__L_c'  131.17      'P'     % I     Isoleucine      ile
         'lys__L_c'  146.19      'P'     % K     Lysine          lys
         'leu__L_c'  131.17      'P'     % L     Leucine         leu
         'met__L_c'  149.21      'P'     % M     Methionine      met
         'asn__L_c'  132.12      'P'     % N     Asparagine      asn
         'pro__L_c'  115.13      'P'     % P     Proline         pro
         'gln__L_c'  146.14      'P'     % Q     Glutamine       gln
         'arg__L_c'  174.2       'P'     % R     Arginine        arg
         'ser__L_c'  105.09      'P'     % S     Serine          ser
         'thr__L_c'  119.12      'P'     % T     Threonine       thr
         'val__L_c'  117.15      'P'     % V     Valine          val
         'trp__L_c'  204.23      'P'     % W     Tryptophan      trp
         'tyr__L_c'  181.19      'P'     % Y     Tyrosine        tyr
                                     };     
%Get main fractions:
[P,X] = getFraction(model,comps,'P',0);
%Add up any remaining components:
bioPos = find(strcmp(model.rxnNames,'biomass pseudoreaction'));
for i = 1:length(model.mets)
    pos = strcmp(comps(:,1),model.mets{i});
    if sum(pos) == 1
        abundance = -model.S(i,bioPos)*comps{pos,2}/1000;
        X         = X + abundance;
    end
end
C = 0;
L = 0;
D = 0;
R = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,X] = getFraction(model,comps,compType,X)
%Define pseudoreaction name:
rxnName = 'protein pseudoreaction';
%Add up fraction:
fractionPos = strcmpi(model.rxnNames,rxnName);
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
X = X + F;
end