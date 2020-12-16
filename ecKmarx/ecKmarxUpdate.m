% ecKmarxUpdate
%
%   Ivan Domenzain, 2020-09-25
%

%Clone the necessary repos:
git('clone https://github.com/SysBioChalmers/GECKO.git')
git('clone https://github.com/SysBioChalmers/Kluyveromyces_marxianus-GEM.git')

%Load kmar model:
model = load('Kluyveromyces_marxianus-GEM/ModelFiles/mat/Kluyveromyces_marxianus-GEM.mat');
model = model.model;
modelVer = model.modelID(strfind(model.modelID,'_v')+1:end);

%Replace scripts in GECKO:
fileNames = dir('scripts');
for i = 1:length(fileNames)
    fileName = fileNames(i).name;
    if ~startsWith(fileName,'.')
        fullName   = ['scripts/' fileName];
        GECKO_path = dir(['GECKO/**/' fileName]);
        GECKO_path = GECKO_path.folder;
        copyfile(fullName,GECKO_path)
    end
end

%Remove unecessary files 
delete('GECKO/databases/prot_abundance.txt') 

%Run GECKO pipeline:
cd GECKO/geckomat
GECKOver = git('describe --tags');
cd get_enzyme_data
updateDatabases
cd ..
[ecModel,ecModel_batch] = enhanceGEM(model,'COBRA','ecKmarx',modelVer(2:end));
cd ../..

%Move model files:
rmdir('model', 's')
movefile GECKO/models/ecKmarx model
save('model/ecKmarx.mat','ecModel')
save('model/ecKmarx_batch.mat','ecModel_batch')

%Save associated versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['GECKO\t' GECKOver '\n']);
fprintf(fid,['Kmarx\t' modelVer '\n']);
fclose(fid);

%Remove the cloned repos:
rmdir('GECKO', 's')
rmdir('Kluyveromyces_marxianus-GEM', 's')