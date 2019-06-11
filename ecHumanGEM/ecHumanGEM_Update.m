% ecHumanGEM_update
%
%   Ivan Domenzain, 2019-06-10
%

%Clone the necessary repos:
git('clone https://github.com/SysBioChalmers/GECKO.git')
git('clone https://github.com/SysBioChalmers/Human-GEM.git')
mkdir model
%Load kmar model:
model = load('Human-GEM/ModelFiles/mat/humanGEM.mat');
model = model.ihuman;
%Replace scripts in GECKO:
replaceFiles('scripts','GECKO/**/');
%Replace databases in GECKO:
replaceFiles('databases','GECKO/databases/');
%Remove unecessary files 
delete('GECKO/databases/prot_abundance.txt') 
%Run GECKO pipeline:
cd GECKO/geckomat
GECKOver = git('describe --tags');
[ecModel,ecModel_batch] = enhanceGEM(model,'ecHumanGEM');
cd ../..

%Move model files:
rmdir('model', 's')
movefile GECKO/models/ecHumanGEM model
save('model/ecHumanGEM.mat','ecModel')
save('model/ecHumanGEM.mat','ecModel_batch')

%Save associated versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['GECKO\t' GECKOver '\n']);
fclose(fid);

%Remove the cloned repos:
 rmdir('GECKO', 's')
 rmdir('Human-GEM', 's')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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