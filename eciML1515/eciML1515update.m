% ecYeastUpdate
%
%   Benjamin J. Sanchez, 2018-10-08
%

%Clone the necessary repos:
git('clone https://github.com/SysBioChalmers/GECKO.git')
git('clone https://github.com/SysBioChalmers/yeast-GEM.git')

%Load yeast model:
cd yeast-GEM
model    = load('ModelFiles/mat/yeastGEM.mat');
model    = model.model;
yeastVer = model.modelID(strfind(model.modelID,'_v')+1:end);
cd ..

%Replace scripts in GECKO:
fileNames = dir('scripts');
for i = 1:length(fileNames)
    fileName = fileNames(i).name;
    if ~strcmp(fileName,'.') && ~strcmp(fileName,'..')
        fullName   = ['scripts/' fileName];
        GECKO_path = dir(['GECKO/**/' fileName]);
        GECKO_path = GECKO_path.folder;
        copyfile(fullName,GECKO_path)
    end
end

%Run GECKO pipeline:
cd GECKO
GECKOver = git('describe --tags');
cd geckomat/get_enzyme_data
updateDatabases;
cd ..
[ecModel,ecModel_batch] = enhanceGEM(model,'COBRA');
cd ../..

%Move model files:
rmdir('model', 's')
movefile GECKO/models/ecYeastGEM model
save('model/ecYeastGEM.mat','ecModel')
save('model/ecYeastGEM_batch.mat','ecModel_batch')

%Save associated versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['GECKO\t' GECKOver '\n']);
fprintf(fid,['yeast-GEM\t' yeastVer '\n']);
fclose(fid);

%Remove the cloned repos:
rmdir('GECKO', 's')
rmdir('yeast-GEM', 's')
