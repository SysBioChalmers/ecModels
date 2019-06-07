% eciML1515Update
%
%  Ivan Domenzain.  2019-06-07
%

%Clone the necessary repos:
% git('clone https://github.com/SysBioChalmers/GECKO.git')
% %Load iML1515 model:
% cd model
% model    = importModel('iML1515.xml');
% cd ..
%Replace scripts in GECKO:
replaceFiles('scripts','GECKO/**/');
%Replace databases in GECKO:
replaceFiles('databases','GECKO/databases/');
%Run GECKO pipeline:
cd GECKO/geckomat
GECKOver = git('describe --tags');
%[ecModel,ecModel_batch,version] = enhanceGEM(model,'COBRA');
cd ../..
%Move model files:
moveModelFiles(name)
save('model/eciML1515.mat','ecModel')
save('model/eciML1515_batch.mat','ecModel_batch')
%Save associated versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['GECKO\t' GECKOver '\n']);
fprintf(fid,['iML1515\t' version '\n']);
fclose(fid);
%Remove the cloned repos:
rmdir('GECKO', 's')
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
        source      = ['GECKO/models/eciML1515/' fileName];
        destination = ['model/' fileName];
        movefile (source,destination);
    end
end
end