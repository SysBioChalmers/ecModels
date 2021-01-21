import tools
import subprocess as sp
from sys import exit
import logging
import os


logging.basicConfig(level=logging.DEBUG)
l = logging.getLogger(__name__)
system = tools.GECKO_VM()


def matlab_command(gem):
    cmd = """
        cd {}geckomat
        model = load('{}');
        model = model{};
        modelname = '{}';
        [ecModel, ecModel_batch] = enhanceGEM(model,'COBRA', modelname, '{}');
        cd ../models
        save([modelname '/' modelname '.mat']','ecModel');
        save([modelname '/' modelname '_batch.mat'], 'ecModel_batch');
        quit
        """.format(system.install_dir('GECKO'), system.mat_file_location(gem), system.mat_model(gem), gem, system.version(gem))
    l.info(cmd)
    output = sp.check_output(['/usr/local/bin/matlab', '-nodisplay -nosplash -nodesktop -r', '"{}"'.format(cmd)])
    return output.decode('utf-8')

def setup_and_run_GECKO(gem):
    system.git_clone('GECKO', system.config['GECKO']['branch'])
    l.info('Merge scripts folder if it exists')
    cmd = """
        fileNames = dir('{}');
        for i = 1:length(fileNames)
            fileName = fileNames(i).name;
            if ~strncmp(fileName, '.', 1) 
                fullName   = ['{}/' fileName];
                GECKO_path = dir(['{}**/' fileName]);
                GECKO_path = GECKO_path.folder;
                copyfile(fullName,GECKO_path)
            end
        end
    """.format(system.scripts(gem), system.scripts(gem), system.install_dir('GECKO'))
    l.info(cmd)
    output = sp.check_output(['/usr/local/bin/matlab', '-nodisplay -nosplash -nodesktop -r', '"{}"'.format(cmd)])
    l.info(output.decode('utf-8'))
    system.cleanup('GECKO', 'models/' + gem)
    sp.check_call(['mkdir', system.install_dir('GECKO') + 'models/' + gem])

    system.git_checkout(gem)
    l.info('Running MATLAB command')
    matlab_output = matlab_command(gem)
    l.info(matlab_output)
    l.info('Copying resulting model files from the GECKO output folder into the current repository')
    sp.check_call(['cp', '-Rf', system.install_dir('GECKO') + 'models/' + gem + '/.', system.JENKINS_WORKSPACE + gem + '/model/'])
    system.git_add_and_pr(gem, matlab_output)
    system.cleanup('GECKO')


l.info('It has begun')
l.info('Checking all dependencies listed in config.ini of type not gem')
system.check_dependencies()

# Run GECKO wherever needed
for gem in system.gems():
    system.cleanup(gem)
    new_version = system.download(gem)
    old_version = system.version(gem)
    if system.HAS_CHANGES or git_version != old_version:
        if system.HAS_CHANGES:
            l.warning('System config has changed, have to run GECKO on {} {}'.format(gem, new_version))
        else:
            l.warning('{} changed from {} to {}'.format(gem, old_version, new_version))

        l.info('Going to run GECKO on {}, saving config file before running'.format(gem))
        system.version(gem, new_version)
        system.save_config()

        setup_and_run_GECKO(gem)

        l.info('Reverting changes on config file made for {}, proceeding to next gem'.format(gem))
        system.version(gem, old_version)
        system.save_config()
    system.cleanup(gem)
