from tools import GECKO_VM
import subprocess as sp
from sys import exit
import logging
import os


logging.basicConfig(level=logging.DEBUG)
l = logging.getLogger(__name__)
system = GECKO_VM()


def matlab_command(gem):
    # Temporary fix, use the devel branch of GECKO
    sp.check_call(['git', 'checkout', '-b', 'devel'], cwd=(system.install_dir('GECKO')))
    cmd = """
        cd geckomat
        model = load('{}');
        model = model.model;
        modelname = '{}';
        [ecModel, ecModel_batch] = enhanceGEM(model,'COBRA', modelname, '{}');
        baseFolder = '{}'
        save([baseFolder modelname '/model/' modelname '.mat']','ecModel');
        save([baseFolder modelname '/model/' modelname '_batch.mat'], 'ecModel_batch');
        quit
        """.format(system.mat_file_location(gem), gem, system.version(gem), system.JENKINS_WORKSPACE)
    print(cmd)
    output = sp.check_call(['/usr/local/bin/matlab', '-nodisplay -nosplash -nodesktop -r', '"{}"'.format(cmd)], cwd=(system.install_dir('GECKO')))
    return output.decode('utf-8')

def setup_and_run_GECKO(gem):
    system.git_clone('GECKO')
    if os.path.exists(system.scripts(gem)) and os.path.exists(system.databases(gem)):
        # Merge scripts folder if it exists
        sp.check_call(['cp', '-Rf', system.scripts(gem), system.install_dir('GECKO') + 'scripts'])
        # Remove the databases in GECKO
        # system.cleanup('GECKO', 'databases')
        # Merge databases folder if it exists
        sp.check_call(['cp', '-Rf', system.databases(gem), system.install_dir('GECKO') + 'databases'])
        # Rm the currently stored ecYeastGEM in the models directory
        # system.cleanup('GECKO', 'models')
    else:
        l.critical('Expected folders for {} are missing, check:\n\t{}\n\t{}'.format(gem, system.scripts(gem), system.databases(gem)))
        return

    system.git_checkout(gem)
    l.info('Running MATLAB command')
    l.info(matlab_command(gem))

    l.info('Copying resulting model files from the GECKO output folder into the current repository')
    sp.check_call(['cp', '-Rf', system.install_dir('GECKO') + 'model/ .'])
    # TODO
    system.git_add_and_pr(gem)
    system.cleanup('GECKO')


l.info('It has begun')
l.info('Checking all dependencies listed in config.ini of type not gem')
system.check_dependencies()

# Run GECKO wherever needed
for gem in system.gems():
    system.cleanup(gem)
    # TODO clone or download
    system.git_clone(gem)
    old_version = system.version(gem)
    git_version = system.git_tag(gem)
    if system.HAS_CHANGES or git_version != old_version:
        if system.HAS_CHANGES:
            l.warning('System config has changed, have to run GECKO on {} {}'.format(gem, git_version))
        else:
            l.warning('{} changed from {} to {}'.format(gem, old_version, git_version))

        l.info('Going to run GECKO on {}, saving config file before running'.format(gem))
        system.version(gem, git_version)
        system.save_config()

        setup_and_run_GECKO(gem)

        l.info('Reverting changes on config file made for {}, proceeding to next gem'.format(gem))
        system.version(gem, old_version)
        system.save_config()
    system.cleanup(gem)
