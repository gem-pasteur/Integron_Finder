#!/usr/bin/env python
# -*- coding: utf-8 -*-


import time
import sys
import os

from ConfigParser import SafeConfigParser, NoSectionError, NoOptionError
try:
    import warnings
except ImportError:
    warnings = None
from distutils import log, dir_util
from distutils.core import setup
from distutils.core import Command
from distutils.dist import Distribution
from distutils.command.build import build
from distutils.command.install import install
from distutils.command.install_data import install_data as _install_data
#from distutils.command.build_scripts import build_scripts as _build_scripts
from distutils.command.install_scripts import install_scripts as _install_scripts
#from distutils.command.sdist import sdist
from distutils.errors import DistutilsFileError, DistutilsOptionError, DistutilsPlatformError
from distutils.versionpredicate import VersionPredicate
from distutils.util import subst_vars as distutils_subst_vars
from distutils.util import get_platform, change_root, convert_path

class check_and_build( build ):

    def run(self):
        chk = True
        for req in require_python:
            chk &= self.check_python(req)
        for req in require_packages:
            chk &= self.check_package(req)
        if not chk:
            sys.exit(1)
        build.run(self)

    def check_python(self, req):
        chk = VersionPredicate(req)
        ver = '.'.join([str(v) for v in sys.version_info[:2]])
        if not chk.satisfied_by(ver):
            log.error("Invalid python version, expected {0}".format(req))
            return False
        return True

    def check_package(self, req):
        chk = VersionPredicate(req)
        try:
            mod = __import__(chk.name)
        except:
            log.error("Missing mandatory {0} python module".format(chk.name))
            return False
        for v in [ '__version__', 'version' ]:
            ver = getattr(mod, v, None)
            break
        try:
            if ver and not chk.satisfied_by(ver):
                log.error("Invalid module version, expected {0}".format(req))
                return False
        except:
            pass
        return True




class install_integron_finder(install):

    #I use record to store all installed files and reuse this record file for uninstall
    #so this option is not available anymore for the users
    for i, opt in enumerate(install.user_options):
        if opt[0] == 'record=':
            install.user_options.pop(i)

    def initialize_options(self):
        install.initialize_options(self)

    def finalize_options(self):
        install.finalize_options(self)
        self.record = self.distribution.uninstall_files
        with open(self.distribution.uninstall_prefix, "w") as _f:
            _f.write("[install]\n")
            _f.write('install_lib = {}\n'.format(os.path.normpath(self.install_lib)))

    def run(self):
        inst = self.distribution.command_options.get('install')
        install.run(self)



class install_scripts(_install_scripts):

    def initialize_options(self):
        _install_scripts.initialize_options(self)

    def finalize_options(self):
        _install_scripts.finalize_options(self)
        with open(self.distribution.uninstall_prefix, "a") as _f:
            _f.write('install_scripts = {}\n'.format(os.path.normpath(self.install_dir)))

    def run(self):
        inst = self.distribution.command_options.get('install')
        vars_2_subst = {'PREFIX': inst['prefix'][1] if 'prefix' in inst else '',
                        #'PREFIXCONF' : os.path.join(get_install_conf_dir(inst), 'integron_finder'),
                        'PREFIXDATA' : os.path.join(get_install_data_dir(inst), 'integron_finder'),
                        #'PREFIXDOC' : os.path.join(get_install_doc_dir(inst), 'integron_finder'),
                        'VERSION' : self.distribution.get_version(),
                        }
        for _file in self.distribution.fix_scripts:
            input_file = os.path.join(self.build_dir, _file)
            output_file = input_file + '.tmp'
            subst_vars(input_file, output_file, vars_2_subst)
            os.unlink(input_file)
            self.move_file(output_file, input_file)
        _install_scripts.run(self)


class install_data(_install_data):

    user_options = [
        ('install-dir=', 'd',
         "base directory for installing data files "
         "(default: installation base dir)"),
        ('root=', None,
         "install everything relative to this alternate root directory"),
        ('force', 'f', "force installation (overwrite existing files)"),
        ]

    boolean_options = ['force']

    def finalize_options(self):
        inst = self.distribution.command_options.get('install')
        self.install_dir = get_install_data_dir(inst)
        self.set_undefined_options('install',
                                   ('root', 'root'),
                                   ('force', 'force'),
                                  )
        self.prefix_data = self.install_dir
        self.files_2_install = self.distribution.data_files
        with open(self.distribution.uninstall_prefix, "a") as _f:
            _f.write('install_data = {}\n'.format(self.install_dir))


    def run(self):
        self.mkpath(self.install_dir)
        for f in self.files_2_install:
            if isinstance(f, str):
                if not os.path.exists(f):
                    log.warn("WARNING the document {} cannot be found, installation skipped".format(f))
                # it's a simple file, so copy it
                f = convert_path(f)
                if self.warn_dir:
                    self.warn("setup script did not provide a directory for "
                              "'{0}' -- installing right in '{1}'".format((f, self.install_dir)))
                (out, _) = self.copy_file(f, self.install_dir)
                self.outfiles.append(out)
            else:
                # it's a tuple with path to install to and a list of path
                dir = convert_path(f[0])
                if not os.path.isabs(dir):
                    dir = os.path.join(self.install_dir, dir)
                elif self.root:
                    dir = change_root(self.root, dir)
                self.mkpath(dir)
                if f[1] == []:
                    # If there are no files listed, the user must be
                    # trying to create an empty directory, so add the
                    # directory to the list of output files.
                    self.outfiles.append(dir)
                else:
                    # Copy files, adding them to the list of output files.
                    for data in f[1]:
                        data = convert_path(data)#return name that will work on the native filesystem
                        if not os.path.exists(data):
                            log.warn("WARNING the document {} cannot be found, installation skipped".format(data))
                            continue
                        if os.path.isdir(data):
                            out = self.copy_tree(data, dir)
                            self.outfiles.extend(out)
                        else:
                            (out, _) = self.copy_file(data, dir)
                            self.outfiles.append(out)



class Uninstall(Command):

    description = "remove installed files"

    user_options = []

    def initialize_options (self):
        self.install_scripts = None
        self.install_lib = None
        self.install_data = None
        #self.install_doc = None


    def finalize_options(self):
        self.parser = SafeConfigParser()
        if not os.path.exists(self.distribution.uninstall_prefix):
            raise DistutilsFileError("Cannot unistall integron_finder.\n{}: No such file".format(self.distribution.uninstall_prefix))
        used_files = self.parser.read(self.distribution.uninstall_prefix)
        for attr in [attr for attr in vars(self) if attr.startswith('install_')]:
            try:
                value = self.parser.get('install', attr)
            except(NoSectionError, NoOptionError):
                continue
            setattr(self, attr, value)

    def run(self):
        prefixes = []
        for attr in [attr for attr in vars(self) if attr.startswith('install_')]:
            prefixes.append( getattr(self, attr))
        print "prefixes = ", prefixes

        def clean_tree(_dir):
            find_prefix = False
            for prefix in prefixes:
                #print "== ", prefix, ".find(",_dir,") = ",prefix.find(_dir)
                if prefix.find(_dir) != -1:
                    find_prefix = True
                    return prefix

            #print "find_prefix =",prefix
            if find_prefix:
                return
            try:
                if not self.dry_run:
                    os.rmdir(_dir)
                log.info("remove dir {}".format(_dir))
            except OSError as err:
                if err.errno == os.errno.ENOTEMPTY:
                    return
                else:
                    self.warn(err)
                    return
            clean_tree(os.path.dirname(_dir))

        try:
            with open(self.distribution.uninstall_files) as record_file:
                for path in record_file:
                    path = os.path.normpath(path.strip())
                    try:
                        if not self.dry_run:
                            os.unlink(path)
                        log.info("remove file {}".format(path))
                    except Exception as err:
                        pass
                    _dir = os.path.dirname(path)
                    clean_tree(_dir)
        except IOError as err:
            msg = "Cannot unistall integron.\n"
            if err.errno == os.errno.ENOENT:
                msg += "Cannot access \"{}\": No such file".format(self.distribution.uninstall_files)
            elif err.errno == os.errno.EACCES:
                msg += "Cannot access \"{}\": Permission denied".format(self.distribution.uninstall_files)
            else:
                msg += str(err)
            raise DistutilsFileError(msg)


class UsageDistribution(Distribution):

    def __init__(self, attrs = None):
        #It's important to define opotions before to call __init__
        #otherwise AttributeError: UsageDistribution instance has no attribute 'conf_files'
        #self.conf_files = None
        #self.doc_files = None
        self.fix_prefix = None
        self.fix_scripts = None
        #self.fix_conf = None
        self.uninstall_prefix = os.path.join(os.path.dirname(__file__), "uninstall.cfg")
        self.uninstall_files = os.path.join(os.path.dirname(__file__), "uninstall_files")
        Distribution.__init__(self, attrs=attrs)
        self.common_usage = """\
Common commands: (see '--help-commands' for more)

  setup.py build      will build the package underneath 'build/'
  setup.py install    will install the package
  setup.py uninstall  will uninstall every installed files
"""


def get_install_data_dir(inst):

    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])

    if 'install_data' in inst:
        install_dir = inst['install_data'][1]
    elif 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'share')
    else:
        install_dir = os.path.join('/', 'usr', 'share')
    return install_dir



def subst_vars(src, dst, vars):
    try:
        src_file = open(src, "r")
    except os.error as err:
        raise DistutilsFileError, "could not open '{0}': {1)".format(src, err)
    try:
        dest_file = open(dst, "w")
    except os.error as err:
        raise DistutilsFileError, "could not create '{0}': {1}".format(dst, err)
    with src_file:
        with dest_file:
            for line in src_file:
                new_line = distutils_subst_vars(line, vars)
                dest_file.write(new_line)



###################################################
#                                                 #
# the configuration of the installer start bellow #
#                                                 #
###################################################

require_python = ['python (>=2.7, <3.0)']

require_packages = ['numpy (>=1.7.0)',
                    'pandas (>=0.18.0)',
                    'Bio (>=1.65)',
                    'matplotlib (>=1.4.2)']

setup(name='integron_finder',
      version="1.4",
      description="""Integron Finder aims at detecting integrons in DNA sequences
       by finding particular features of the integron:
       the attC sites, the integrase
       and when possible attI site and promoters.""",
      author="Jean Cury",
      author_email="jean.cury@normalesup.org",
      classifiers=[
                     'Operating System :: POSIX',
                     'Programming Language :: Python',
                     'Topic :: Bioinformatics' ,
                    ],
      scripts=['integron_finder'],
      #(dataprefix +'where to put the data in the install, [where to find the data in the tar ball]
      data_files=[('integron_finder/Functional_annotation', ['data/Functional_annotation/']),
                  ('integron_finder/Models', ['data/Models/']),
                  ],
      #file where some variable must be fix by install_conf
      #fix_conf=['etc/integron_finder.conf'],
      #file where some variable must be fix by integron_finder_install
      #fix_prefix=['integron_finder'],
      fix_scripts=['integron_finder'],
      cmdclass={'build' : check_and_build,
                'install' : install_integron_finder,
                'install_scripts' : install_scripts,
                'install_data' : install_data,
                'uninstall' : Uninstall,
              },
      distclass=UsageDistribution
      )
