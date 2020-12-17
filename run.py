from tools import GECKO_VM
import subprocess as sp
from sys import exit
import logging
import os


CONFIG_HAS_UPDATES = False

logging.basicConfig(level=logging.DEBUG)
l = logging.getLogger(__name__)
system = GECKO_VM()


def run_GECKO(gem):
    system.cleanup('GECKO', 'models')
    l.critical('running matlab here')
    # also copy the output models from gecko to the right gem folder

def setup_and_run_GECKO(gem):
    system.git_clone('GECKO')
    # Make sure that gem output dir exists
    output_dir = system.output_dir(gem)
    if not os.path.exists(output_dir):
        l.critical('Trying to run GECKO on gem {} provided by config, but required folder does not exist'.format(gem))
        go to next gem if this crashes
    # Merge scripts and databases if they exist
    if os.path.exists(system.scripts(gem)):
        go to next gem if this crashes
        sp.check_call(['cp', '-Rf', system.scripts(gem), system.install_dir('GECKO') + 'scripts'])
    if os.path.exists(system.databases(gem)):
        go to next gem if this crashes
        sp.check_call(['cp', '-Rf', system.databases(gem), system.install_dir('GECKO') + 'databases'])
    system.git_checkout(gem)
    run_GECKO(gem)
    system.git_add_and_pr(gem)
    system.cleanup('GECKO')


l.info("It has begun")

# Check Matlab version
# cmd = sp.check_output(['matlab', '-nodisplay -nojvm -nosplash -nodesktop -r', '"disp(version); quit"'])
# m_version = ' '.join(cmd.decode('utf-8').split()[-3:-1])
# if m_version != system.version('MATLAB'):
#     CONFIG_HAS_UPDATES = True
#     l.warning('MATLAB changed from {} to {}'.format(system.version('MATLAB'), m_version))
#     system.version('MATLAB', m_version)

# Check libSBML version
with open(system.install_dir('libSBML') + 'VERSION.txt') as f:
    l_version = f.readline().strip()
    if system.version('libSBML') != l_version:
        CONFIG_HAS_UPDATES = True
        l.warning('libSMBL changed from {} to {}'.format(system.version('libSBML'), l_version))
        system.version('libSBML', l_version)

# Check COBRA, RAVEN, GECKO versions
# system.cleanup('GECKO')
# system.git_clone('GECKO')
# for tool in ['COBRA', 'RAVEN', 'GECKO']:
#     tool_version = system.git_tag(tool)
#     if tool_version != system.version(tool):
#         l.warning('{} changed from {} to {}'.format(tool, system.version(tool), tool_version))
#         CONFIG_HAS_UPDATES = True
#         system.version(tool, tool_version)
# Cleanup dummy GECKO install
system.cleanup('GECKO')

# Run GECKO wherever needed
for gem in system.gems():
    system.cleanup(gem)
    system.git_clone(gem)
    g_version = system.git_tag(gem)
    if CONFIG_HAS_UPDATES or g_version != system.version(gem):
        l.warning('{} changed from {} to {}'.format(gem, system.version(gem), g_version))
        system.version(gem, g_version)
        setup_and_run_GECKO(gem)
        # revert version to not affect PR of the next gem
    system.cleanup(gem)
