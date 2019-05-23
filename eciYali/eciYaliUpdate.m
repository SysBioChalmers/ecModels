% eciYaliUpdate
%
%   Ivan Domenzain, 2019-02-11
%

%Clone the necessary repos:
git('clone https://github.com/SysBioChalmers/GECKO.git')
git('clone https://github.com/SysBioChalmers/Yarrowia_lipolytica_W29-GEM.git')

%Load kmar model:
model    = load('Yarrowia_lipolytica_W29-GEM/ModelFiles/mat/iYali.mat');
model    = model.model;
modelVer = model.description(strfind(model.description,'_v')+1:end);

%Replace scripts in GECKO:
replaceFiles('scripts','GECKO/**/');
%Replace databases in GECKO:
replaceFiles('databases','GECKO/databases/');
%Remove unecessary files
delete('GECKO/databases/prot_abundance.txt')
%Run GECKO pipeline:
cd GECKO/geckomat
GECKOver = git('describe --tags');
[ecModel,ecModel_batch] = enhanceGEM(model,'COBRA','eciYali');
cd ../..

%Move model files:
rmdir('model', 's')
movefile GECKO/models/eciYali model
save('model/eciYali.mat','ecModel')
save('model/eciYali_batch.mat','ecModel_batch')

%Save associated versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['GECKO\t' GECKOver '\n']);
fprintf(fid,['iYali\t' modelVer '\n']);
fclose(fid);

%Remove the cloned repos:
 rmdir('GECKO', 's')
 rmdir('Yarrowia_lipolytica_W29-GEM', 's')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function replaceFiles(fileType,path)
fileNames = dir(fileType);
for i = 1:length(fileNames)
    fileName = fileNames(i).name;
    if ~strcmp(fileName,'.') && ~strcmp(fileName,'..') && ~strcmp(fileName,'.DS_Store')
        fullName   = [fileType '/' fileName];
        GECKO_path = dir([path fileName]);
        GECKO_path = GECKO_path.folder;
        copyfile(fullName,GECKO_path)
    end
end
end