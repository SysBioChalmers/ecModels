from configparser import ConfigParser
import subprocess as sp
from sys import exit
import errno
import logging
from os import getcwd, environ
import time
import datetime


# Constants
CONFIGFILE = 'config.ini'
URL = 'url'
IDIR = 'install_dir'
SCRIPTSDIR = 'scripts'
DBSDIR = 'databases'

l = logging.getLogger(__name__)


class GECKO_VM:
    """Interact with configuration variables."""

    HAS_CHANGES = False
    config = ConfigParser()
    JENKINS_WORKSPACE = getcwd() + "/"

    def __init__(self):
        self.config.read(CONFIGFILE)
        if len(self.config.sections()) == 0:
            l.critical('Missing files, check {}'.format(CONFIGFILE))
            exit(errno.ENOENT)

    def version(self, section, new_version=None):
        if new_version:
            self.config[section]['version'] = new_version
        else:
            return self.config[section]['version']

    def base_dir(self):
        return self.config['BASE'][IDIR]

    def install_dir(self, section):
        return self.config['BASE'][IDIR] + self.config[section][IDIR]

    def gems(self):
        gems = []
        for section in self.config.sections():
            if self.config.has_option(section, 'type'):
                gems.append(section)
        return gems

    def cleanup(self, section, subdir=''):
        sp.check_call(['rm', '-rf', self.install_dir(section) + subdir])
        l.info('Have removed {}'.format(self.install_dir(section) + subdir))

    def scripts(self, gem):
        return '{}/{}'.format(gem, SCRIPTSDIR)

    def databases(self, gem):
        return '{}/{}'.format(gem, DBSDIR)

    def mat_file_location(self, gem):
        return self.install_dir(gem) + "/ModelFiles/mat/" + self.config[gem]['mat_filename']

    def pr_target(self):
        return self.config['BASE']['pull_request_target']

    def git_clone(self, section):
        sp.check_call(['git', 'clone', '--depth=1', self.config[section][URL], self.install_dir(section)], stdout=sp.DEVNULL, stderr=sp.STDOUT)

    def git_tag(self, thing):
        cmd = sp.Popen(['git', 'describe', '--tags'], cwd=self.install_dir(thing), stdout=sp.PIPE)
        tool_version, _ = cmd.communicate()
        return tool_version.decode('utf-8').strip()

    def git_checkout(self, gem):
        sp.check_call(['git', 'checkout', '-B', self.__branch_name(gem)])

    def git_add_and_pr(self, gem, matlab_output):
        sp.check_call(['git', 'add', gem])
        sp.check_call(['git', 'add', 'config.ini'])
        try:
            # If nothing was addded (no changes) the commit will exit with an error so we can delete the branch
            sp.check_call(['git', 'commit', '-m', '"chore: update {} based on {}"'.format(gem, self.version(gem))], stdout=sp.DEVNULL, stderr=sp.STDOUT)
            l.critical('Will push and create PR')
            # Create PR and also push
            pr_filename = "/tmp/githubpr"
            with open(pr_filename, "w") as f:
                f.write("update {} based on {}\n\n".format(gem, self.version(gem)))
                f.write(matlab_output)
            my_env = environ.copy()
            sp.check_call(['hub', 'pull-request', '--file', pr_filename, '-b', self.pr_target(), '-p'], env=my_env)
        except sp.CalledProcessError:
            l.warning('While upgrading {} to {} no changes were detected, checking out {}'.format(gem, self.version(gem), self.pr_target()))
            sp.check_call(['git', 'checkout', '-f', self.pr_target()], stdout=sp.DEVNULL, stderr=sp.STDOUT)

    def check_dependencies(self):
        # Check Matlab version
        cmd = sp.check_output(['/usr/local/bin/matlab', '-nodisplay -nosplash -nodesktop -r', '"disp(version); quit"'])
        m_version = ' '.join(cmd.decode('utf-8').split()[-3:-1])
        if m_version != self.version('MATLAB'):
            self.HAS_CHANGES = True
            l.warning('MATLAB changed from {} to {}'.format(self.version('MATLAB'), m_version))
            self.version('MATLAB', m_version)
        else:
            l.info('MATLAB is still {}'.format(m_version))
        # Check libSBML, Gurobi version; these version files have been manually created
        for tool in ['libSBML', 'Gurobi']:
            with open(self.install_dir(tool) + 'VERSION.txt') as f:
                tool_version = f.readline().strip()
                if self.version(tool) != tool_version:
                    self.HAS_CHANGES = True
                    l.warning('{} changed from {} to {}'.format(tool, self.version(tool), tool_version))
                    self.version(tool, tool_version)
                else:
                    l.info('{} is still {}'.format(tool, tool_version))
        # Check COBRA, RAVEN, GECKO versions
        self.cleanup('GECKO')
        self.git_clone('GECKO')
        for tool in ['COBRA', 'RAVEN', 'GECKO']:
            tool_version = self.git_tag(tool)
            if tool_version != self.version(tool):
                l.warning('{} changed from {} to {}'.format(tool, self.version(tool), tool_version))
                self.HAS_CHANGES = True
                self.version(tool, tool_version)
            else:
                l.info('{} is still {}'.format(tool, tool_version))
        # Cleanup dummy GECKO install
        self.cleanup('GECKO')

    def save_config(self):
        with open(CONFIGFILE, 'w') as configfile:
            self.config.write(configfile)

    def __branch_name(self, gem):
        time_seconds = time.time()
        timestamp = datetime.datetime.fromtimestamp(time_seconds).strftime('%Y-%m-%d-%H-%M-%S')
        return 'update/{}/{}/{}'.format(gem, self.version(gem), timestamp)
