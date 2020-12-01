% ecHumanGEM_update
%
%   Ivan Domenzain, 2020-12-01
%

%Clone the necessary repos:
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
git('checkout devel')
git('pull origin')
cd ..
git('clone https://github.com/SysBioChalmers/Human-GEM.git')
%Load human-based model:
load('Human-GEM/model/Human-GEM.mat')
model     = ihuman;
modelName = 'ecHumanGEM';

if isfield(model,'version')
    modelVer = model.version;
end
%Replace scripts in GECKO:
fileNames = dir('scripts');
for i = 1:length(fileNames)
    fileName = fileNames(i).name;
    if ~strcmp(fileName,'.') && ~strcmp(fileName,'..') && ~strcmp(fileName,'.DS_Store')
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
[ecModel,ecModel_batch] = enhanceGEM(model,'COBRA',modelName,modelVer);
cd ../..
%Move model files:
rmdir('model', 's')
movefile(['GECKO/models/' modelName],'model')
save(['model/' modelName '.mat'],'ecModel')
save(['model/' modelName '_batch.mat'],'ecModel_batch')

%Save associated versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['GECKO\t' GECKOver '\n']);
fprintf(fid,[modelName '\t' modelVer '\n']);
fclose(fid);

%Remove the cloned repos:
 rmdir('GECKO', 's')
 rmdir('Human-GEM', 's')