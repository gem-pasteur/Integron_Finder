#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
if sys.version_info[0] == 3:
    sys.exit("Sorry, Python 3 is not supported yet")

import os

try:
    import warnings
except ImportError:
    warnings = None

from distutils import log
from distutils.command.install_data import install_data as _install_data
from distutils.errors import DistutilsFileError
from distutils.util import subst_vars as distutils_subst_vars
from distutils.util import change_root, convert_path

from setuptools import setup
from setuptools.dist import Distribution
from setuptools.command.install import install as _install
from setuptools.command.install_scripts import install_scripts as _install_scripts


class install_scripts(_install_scripts):

    def finalize_options(self):
        inst = self.distribution.command_options.get('install')
        inst = {} if inst is None else inst
        _install_scripts.finalize_options(self)

    def run(self):
        if os.path.exists(self.build_dir):
            installer = 'pip'
        else:
            installer = 'fucking_setuptools'

        script_tmp_dir = self.build_dir if installer != 'fucking_setuptools' else self.install_dir

        def subst_file(_file, vars_2_subst):
            input_file = os.path.join(script_tmp_dir, _file)
            output_file = input_file + '.tmp'
            subst_vars(input_file, output_file, vars_2_subst)
            os.unlink(input_file)
            self.move_file(output_file, input_file)

        # if setup.py install is used without using setup.py build
        # setuptools creates a directory build/bdist.linux-x86_64/egg/EGG-INFO/scripts/
        # but not build/scripts
        # and without attribute to grasp this directory :-(((
        # so we need to run _install_scripts first
        # to create build/scripts (self.build_dir)
        # and then substitute variables in this dir

        if installer == 'fucking_setuptools':
            _install_scripts.run(self)

        inst = self.distribution.command_options.get('install')
        inst = {} if inst is None else inst
        if self.distribution.fix_scripts is not None:
            vars_2_subst = {'PREFIX': inst['prefix'][1] if 'prefix' in inst else '',
                            'PREFIXDATA': os.path.join(get_install_data_dir(inst), 'integron_finder'),
                            'VERSION': self.distribution.get_version(),
                            }
            for _file in self.distribution.fix_scripts:
                subst_file(_file, vars_2_subst)
                pass
        if installer == 'pip':
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
        inst = {} if inst is None else inst
        self.install_dir = get_install_data_dir(inst)
        self.set_undefined_options('install',
                                   ('root', 'root'),
                                   ('force', 'force'),
                                  )
        self.prefix_data = self.install_dir
        self.files_2_install = self.distribution.data_files


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
                              "'{0}' -- installing right in '{1}'".format(f, self.install_dir))
                (out, _) = self.copy_file(f, self.install_dir)
                self.outfiles.append(out)
            else:
                # it's a tuple with path to install to and a list of path
                dir_ = convert_path(f[0])
                if not os.path.isabs(dir_):
                    dir_ = os.path.join(self.install_dir, dir_)
                elif self.root:
                    dir_ = change_root(self.root, dir_)
                self.mkpath(dir_)
                if f[1] == []:
                    # If there are no files listed, the user must be
                    # trying to create an empty directory, so add the
                    # directory to the list of output files.
                    self.outfiles.append(dir_)
                else:
                    # Copy files, adding them to the list of output files.
                    for data in f[1]:
                        data = convert_path(data)  # return name that will work on the native filesystem
                        if not os.path.exists(data):
                            log.warn("WARNING the document {} cannot be found, installation skipped".format(data))
                            continue
                        if os.path.isdir(data):
                            out = self.copy_tree(data, dir_)
                            self.outfiles.extend(out)
                        else:
                            (out, _) = self.copy_file(data, dir_)
                            self.outfiles.append(out)


class install_doc(install_data):

    _install.sub_commands += [('install_doc', lambda self: not self.no_doc)]

    description = "installation directory for documentation files"

    setattr(_install, 'install_doc', None)
    setattr(_install, 'no_doc', None)

    _install.user_options.append(('install-doc=', None, description))
    _install.user_options.append(('no-doc', None, 'do not install documentation'))

    user_options = [
        ('install-doc=', 'd', "base directory for installing documentation files "
                              "(default: installation base dir share/doc)"),
        ('root=', None, "install everything relative to this alternate root directory"),
        ('force', 'f', "force installation (overwrite existing files)"),
        ('no-doc', None, 'do not install documentation')
        ]

    boolean_options = ['force']

    def initialize_options(self):
        install_data.initialize_options(self)
        self.install_doc = None
        self.no_doc = None
        self.files_2_install = self.distribution.doc_files

    def finalize_options(self):
        inst = self.distribution.command_options.get('install')
        inst = {} if inst is None else inst
        self.install_dir = get_install_doc_dir(inst)
        self.set_undefined_options('install',
                                   ('root', 'root'),
                                   ('force', 'force'),
                                  )
        self.prefix_data = self.install_dir
        self.no_doc = inst.get('no_doc', ('command line', False))[1]

    def run(self):
        install_data.run(self)


class UsageDistribution(Distribution):

    def __init__(self, attrs=None):
        #It's important to define potions before to call __init__
        #otherwise AttributeError: UsageDistribution instance has no attribute 'conf_files'
        self.doc_files = None
        self.fix_prefix = None
        self.fix_scripts = None
        Distribution.__init__(self, attrs=attrs)
        self.common_usage = """\
Common commands: (see '--help-commands' for more)

  setup.py build      will build the package underneath 'build/'
  setup.py install    will install the package
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


def get_install_doc_dir(inst):
    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])

    if 'install_doc' in inst:
        install_dir = inst['install_doc'][1]
    elif 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'share', 'integron_finder', 'doc')
    else:
        install_dir = os.path.join('/', 'usr', 'share', 'integron_finder', 'doc')
    return install_dir


def subst_vars(src, dst, vars):
    try:
        src_file = open(src, "r")
    except os.error as err:
        raise DistutilsFileError("could not open '{0}': {1)".format(src, err))
    try:
        dest_file = open(dst, "w")
    except os.error as err:
        raise DistutilsFileError("could not create '{0}': {1}".format(dst, err))
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

setup(name='integron_finder',
      version="1.5",
      description="Integron Finder aims at detecting integrons in DNA sequences",
      long_description="""Integron Finder aims at detecting integrons in DNA sequences
by finding particular features of the integron:
  - the attC sites
  - the integrase
  - and when possible attI site and promoters.""",
      author="Jean Cury",
      author_email="jean.cury@normalesup.org",
      url="https://github.com/gem-pasteur/Integron_Finder/",
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 2',
          'Intended Audience :: Science/Research',
          'Scientific/Engineering :: Bio-Informatics'
          ],
      install_requires=['numpy>=1.7.0',
                        'matplotlib>=1.4.2',
                        'pandas>=0.18.0',
                        'biopython>=1.65',
                        ],

      scripts=['integron_finder'],
      #file where some variable must be fix by install
      fix_scripts=['integron_finder'],

      #(dataprefix +'where to put the data in the install, [where to find the data in the tar ball]
      data_files=[('integron_finder/data/Functional_annotation', ['data/Functional_annotation/']),
                  ('integron_finder/data/Models', ['data/Models/']),
                  ],
      doc_files=[('html', ['doc/build/html/']),
                 ('pdf', ['doc/build/latex/IntegronFinder.pdf'])],
      cmdclass={
                'install_scripts': install_scripts,
                'install_data': install_data,
                'install_doc': install_doc,
              },
      distclass=UsageDistribution
      )


