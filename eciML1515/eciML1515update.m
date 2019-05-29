% eciML1515Update
%
%  Ivan Domenzain.  2019-05-28
%

%Clone the necessary repos:
git('clone https://github.com/SysBioChalmers/GECKO.git')
%Load iML1515 model:
cd model
model    = importModel('iML1515.xml');
cd ..
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

%Run GECKO pipeline:
cd GECKO
GECKOver = git('describe --tags');
cd geckomat
[ecModel,ecModel_batch,version] = enhanceGEM(model,'COBRA');
cd ../..

%Move model files:
rmdir('model', 's')
movefile GECKO/models/eciML1515 model
save('model/eciML1515.mat','ecModel')
save('model/eciML1515_batch.mat','ecModel_batch')

%Save associated versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['GECKO\t' GECKOver '\n']);
fprintf(fid,['iML1515\t' version '\n']);
fclose(fid);

%Remove the cloned repos:
rmdir('GECKO', 's')
