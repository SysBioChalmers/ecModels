% ecKmarxUpdate
%
%   Ivan Domenzain, 2019-02-06
%

%Clone the necessary repos:
git('clone https://github.com/SysBioChalmers/GECKO.git')
git('clone https://github.com/SysBioChalmers/Kluyveromyces_marxianus-GEM.git')

%Load kmar model:
model = load('Kluyveromyces_marxianus-GEM/ModelFiles/mat/kmar.mat');
model = model.model;
%Replace scripts in GECKO:
replaceFiles('scripts','GECKO/**/');
%Replace databases in GECKO:
replaceFiles('databases','GECKO/databases/');
%Remove unecessary files 
delete('GECKO/databases/prot_abundance.txt') 
%Run GECKO pipeline:
cd GECKO/geckomat
GECKOver = git('describe --tags');
cd get_enzyme_data
updateDatabases('kmx')
cd ..
[ecModel,ecModel_batch,version] = enhanceGEM(model,'COBRA','ecKmarx');
cd ../..
%Move model files:
moveModelFiles('ecKmarx')
save('model/ecKmarx.mat','ecModel')
save('model/ecKmarx_batch.mat','ecModel_batch')
%Save associated versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['GECKO\t' GECKOver '\n']);
fprintf(fid,['Kmarx\t' version '\n']);
fclose(fid);
%Remove the cloned repos:
 rmdir('GECKO', 's')
 rmdir('Kluyveromyces_marxianus-GEM', 's')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function moveModelFiles(name)
cd GECKO/models
fileNames = dir(name);
cd ../..
for i=1:length(fileNames)
    fileName = fileNames(i).name;
    if ~strcmp(fileName,'.') && ~strcmp(fileName,'..') && ~strcmp(fileName,'.DS_Store')
        source      = ['GECKO/models/ecKmarx/' fileName];
        destination = ['model/' fileName];
        movefile (source,destination);
    end
end
end