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
%        id         MW [g/mol]  class     name
comps = {'s_0955'	89.09       'P'     % A     Alanine         ala
         's_0981'	121.16      'P'     % C     Cysteine        cys
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
         's_0509'   221.21      'C'     % chitin
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
         's_0666'   396.65      'L'     % Ergosterol
         's_3656'   633.058     'L'     % ergosteryl palmitoleate
         's_3658'   661.112     'L'     % ergosteryl oleate
         's_3710'   385.306     'L'     % phosphatidyl-L-serine
         'm1640'    867.22      'L'     % triglyceride (average in Yeast; %https://genome.cshlp.org/content/suppl/2003/02/03/13.2.244.DC1/3.pdf)
         'm1641'    259.406     'L'     % fatty acid (average in Yeast; %https://genome.cshlp.org/content/suppl/2003/02/03/13.2.244.DC1/3.pdf)
         'm1648'    470.213     'L'     % 1-phosphatidyl-1D-myo-inositol
         'm1651'    1466.086    'L'     % cardiolipin
         'm1705'    226.07      'L'     % phosphatidate
         'm1700'    269.146     'L'     % phosphatidylethanolamine
         'm1700'    312.085     'L'};     % phosphatidylcholine
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
