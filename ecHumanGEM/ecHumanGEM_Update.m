% ecHumanGEM_update
%
%   Ivan Domenzain, 2019-06-16
%

%Clone the necessary repos:
git('clone https://github.com/SysBioChalmers/GECKO.git')
mkdir model
%Load human-based model model:
model     = loadModel(modelName);
modelName = ['ec' modelName];
mkdir (['model/' modelName])
humanVer  = '';
if isfield(model,'version')
    humanVer = model.version;
else
    while isempty(humanVer)
        humanVer = input('Please enter the model version: ','s');
    end
end
cd ..
%Replace scripts in GECKO:
replaceFiles('scripts','GECKO/**/');
%Replace databases in GECKO:
replaceFiles('databases','GECKO/databases/');
%Remove unecessary files 
delete('GECKO/databases/prot_abundance.txt') 
%Run GECKO pipeline:
cd GECKO/geckomat
GECKOver = git('describe --tags');
[ecModel,ecModel_batch] = enhanceGEM(model,modelName,humanVer);
cd ../..

%Move model files:
rmdir('model', 's')
movefile(['GECKO/models/' modelName],['models/' modelName])
save(['models/' modelName '/' modelName '.mat'],'ecModel')
save(['models/' modelName '/' modelName '_batch.mat'],'ecModel_batch')

%Save associated versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['GECKO\t' GECKOver '\n']);
fprintf(fid,['Human-GEM\t' humanVer '\n']);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = loadModel(name)
if strcmpi(name,'HumanGEM')
    git('clone https://github.com/SysBioChalmers/Human-GEM.git')
    model    = load('Human-GEM/ModelFiles/mat/humanGEM.mat');
    model    = model.ihuman;
else
    cd cell_lines
    try 
        model = load([name '.mat']);
        model = model.model;
        if isfield(model,'rxnFrom')
            model = rmfield(model,'id');
        end
    catch
        error('Model not present in cell_lines')
    end
end
end