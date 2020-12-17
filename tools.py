from configparser import ConfigParser
import subprocess as sp
from sys import exit
import errno
import logging


CONFIGFILE = 'config.ini'
URL = 'url'
IDIR = 'install_dir'
ODIR = 'output_dir'
SCRIPTSDIR = 'scripts'
DBSDIR = 'databases'

l = logging.getLogger(__name__)


class GECKO_VM:
    """Interact with configuration variables."""

    config = ConfigParser()

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
        l.info('Have removed {}'.format(section))

    def output_dir(self, gem):
        return self.base_dir() + self.config[gem][ODIR]

    def scripts(self, gem):
        return self.base_dir() + self.config[gem][ODIR] + SCRIPTSDIR

    def databases(self, gem):
        return self.base_dir() + self.config[gem][ODIR] + DBSDIR

    def git_clone(self, section):
        sp.check_call(['git', 'clone', '--depth=1', self.config[section][URL], self.install_dir(section)], stdout=sp.DEVNULL, stderr=sp.STDOUT)

    def git_tag(self, thing):
        cmd = sp.Popen(['git', 'describe', '--tags'], cwd=self.install_dir(thing), stdout=sp.PIPE)
        tool_version, _ = cmd.communicate()
        return tool_version.decode('utf-8').strip()

    def git_checkout(self, gem):
        sp.check_call(['git', 'checkout', '-b', self.__branch_name(gem)], cwd=self.base_dir())

    def git_add_and_pr(self, gem):
        sp.check_call(['git', 'add', self.output_dir(gem)], cwd=self.base_dir())
        sp.check_call(['git', 'add', 'config.ini'], cwd=self.base_dir())
        try:
            # If nothing was addded (no changes) the commit will exit with an error so we can delete the branch
            sp.check_call(['git', 'commit', '-m', '"chore: update {} based on {}"'.format(gem, self.version(gem))], cwd=self.base_dir(), stdout=sp.DEVNULL, stderr=sp.STDOUT)
            l.critical('WILL PUSH AND PR')
            # sp.check_call(['git', 'push', '--set-upstream', 'origin', self.__branch_name(gem)])
            # sp.check_call(['git', 'request-pull', self.__branch_name(gem), self.configp['BASE'][URL], 'master'])
        except sp.CalledProcessError:
            l.warning('While upgrading {} to {} no changes were detected. Checking out master and deleting temporrary branch.'.format(gem, self.version(gem)))
            sp.check_call(['git', 'checkout', 'master'], cwd=self.base_dir(), stdout=sp.DEVNULL, stderr=sp.STDOUT)
            sp.check_call(['git', 'branch', '-D', self.__branch_name(gem)], cwd=self.base_dir(), stdout=sp.DEVNULL, stderr=sp.STDOUT)

    def __branch_name(self, gem):
        return 'update/{}/{}'.format(gem, self.version(gem))
