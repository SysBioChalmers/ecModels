from tools import GECKO_VM
import subprocess as sp
from sys import exit
import logging
import os


logging.basicConfig(level=logging.DEBUG)
l = logging.getLogger(__name__)
system = GECKO_VM()


def matlab_command(gem):
    cmd = """
        model = load({});
        model = model.model;
        modelname = '{}'
        [ecModel, ecModel_batch] = enhanceGEM(model,'COBRA', modelname);
        movefile GECKO/models/ecYeastGEM model
        cd
        save([modelname '.mat'], 'ecModel')
        save([modelname '_batch.mat'], 'ecModel_batch')
        """.format(system.mat_file_location(gem), gem)


def setup_and_run_GECKO(gem):
    system.git_clone('GECKO')
    if os.path.exists(system.scripts(gem)) and os.path.exists(system.databases(gem)):
        # Merge scripts folder if it exists
        sp.check_call(['cp', '-Rf', system.scripts(gem), system.install_dir('GECKO') + 'scripts'])
        # Remove prot_abundance file from the databases in GECKO
        system.cleanup('GECKO', 'databases/prot_abundance.txt')
        # Merge databases folder if it exists
        sp.check_call(['cp', '-Rf', system.databases(gem), system.install_dir('GECKO') + 'databases'])
        # Rm the currently stored ecYeastGEM in the models directory
        system.cleanup('GECKO', 'models')

        system.git_checkout(gem)
        l.critical('running matlab here')
        also copy the output models from gecko to the right gem folder

        system.git_add_and_pr(gem)
    else:
        l.critical('Expected folders for {} are missing, skipping' gem)
        return
    system.cleanup('GECKO')


l.info("It has begun")

# Check all dependencies listed in config.ini of type not gem
system.check_dependencies()

# Run GECKO wherever needed
for gem in system.gems():
    system.cleanup(gem)
    system.git_clone(gem)
    old_version = system.version(gem)
    git_version = system.git_tag(gem)
    if system.HAS_CHANGES or git_version != old_version:
        if system.HAS_CHANGES:
            l.warning('System config has changed, running GECKO on {} {}'.format(gem, git_version))
        else:
            l.warning('{} changed from {} to {}'.format(gem, old_version, git_version))
        # Make sure to save the config for each branch before adding it to the repo
        system.version(gem, git_version)
        self.save_config()

        setup_and_run_GECKO(gem)
        # Revert gem version to not affect PR of the next gem
        system.version(gem, old_version)
        self.save_config()
    system.cleanup(gem)
