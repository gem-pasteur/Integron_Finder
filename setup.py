#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
if sys.version_info[0] == 3:
    sys.exit("Sorry, Python 3 is not supported yet")

import os
import sysconfig

try:
    import warnings
except ImportError:
    warnings = None

from distutils.errors import DistutilsFileError
from distutils.util import subst_vars as distutils_subst_vars

from setuptools import setup
from setuptools.dist import Distribution
from setuptools.command.install_scripts import install_scripts as _install_scripts


class install_scripts(_install_scripts):

    def finalize_options(self):
        inst = self.distribution.command_options.get('install', {})
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

        inst = self.distribution.command_options.get('install', {})
        if self.distribution.fix_scripts is not None:
            vars_2_subst = {'PREFIX': inst['prefix'][1] if 'prefix' in inst else '',
                            'PREFIXSHARE': os.path.join(get_install_data_dir(inst), 'integron_finder'),
                            'VERSION': self.distribution.get_version(),
                            }
            for _file in self.distribution.fix_scripts:
                subst_file(_file, vars_2_subst)
                pass
        if installer == 'pip':
            _install_scripts.run(self)


class UsageDistribution(Distribution):

    def __init__(self, attrs=None):
        #It's important to define potions before to call __init__
        #otherwise AttributeError: UsageDistribution instance has no attribute 'conf_files'
        #self.doc_files = None
        self.fix_prefix = None
        self.fix_scripts = None
        Distribution.__init__(self, attrs=attrs)
        self.common_usage = """\
Common commands: (see '--help-commands' for more)

  setup.py build      will build the package underneath 'build/'
  setup.py install    will install the package
"""


def get_install_data_dir(inst):
    """
    :param inst: installation option
    :type inst: dict
    :return: the prefix where to install data
    :rtype: string
    """

    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])
    elif 'user' in inst:
        import site
        inst['prefix'] = ('command line', site.USER_BASE)

    if 'install_data' in inst:
        install_dir = inst['install_data'][1]
    elif 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'share')
    else:
        install_dir = os.path.join(sysconfig.get_path('data'), 'share')
    return install_dir


def subst_vars(src, dst, vars):
    """
    substitute variables (string starting with $) in file
    :param src: the file containing variable to substitute
    :type src: string
    :param dst: the destination file
    :type dst: string
    :param vars: the variables to substitute in dict key are variable name
    :type vars: dict
    """
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


def expand_data(data_to_expand):
    """
    From data structure like setup.py data_files (see http://)
     [(directory/where/to/copy/the/file, [path/to/file/in/archive/file1, ...]), ...]
    but instead of the original struct this one accept to specify a directory in elements to copy.

    This function will generate one entry for all *content* of the directory and subdirectory
    recursively, to in fine copy the tree in archive in dest on the host

    the first level of directory itself is not include (which allow to rename it)
    :param data_to_expand:
    :type  data_to_expand: list of tuple
    :return: list of tuple
    """
    def remove_prefix(prefix, path):
        prefix = os.path.normpath(prefix)
        path = os.path.normpath(path)
        to_remove = len([i for i in prefix.split(os.path.sep) if i])
        truncated = [i for i in path.split(os.path.sep) if i][to_remove:]
        truncated = os.path.sep.join(truncated)
        return truncated

    data_struct = []
    for base_dest_dir, src in data_to_expand:
        base_dest_dir = os.path.normpath(base_dest_dir)
        for one_src in src:
            if os.path.isdir(one_src):
                for path, _, files in os.walk(one_src):
                    if not files:
                        continue
                    path_2_create = remove_prefix(one_src, path)
                    data_struct.append(
                        (os.path.join(base_dest_dir, path_2_create), [os.path.join(path, f) for f in files]))
            if os.path.isfile(one_src):
                data_struct.append((base_dest_dir, [one_src]))
    return data_struct


###################################################
#                                                 #
# the configuration of the installer start bellow #
#                                                 #
###################################################

setup(name='integron_finder',
      version="1.5.1.RC1",
      description="Integron Finder aims at detecting integrons in DNA sequences",
      long_description="""Integron Finder aims at detecting integrons in DNA sequences
by finding particular features of the integron:
  - the attC sites
  - the integrase
  - and when possible attI site and promoters.

integron_finder has been describe in a paper published in Nucleic Acid Research.

Identification and analysis of integrons and cassette arrays in bacterial genomes.
Jean Cury; Thomas Jove; Marie Touchon; Bertrand Néron; Eduardo PC Rocha
Nucleic Acids Research 2016; doi:10.1093/nar/gkw319
""",
      author="Jean Cury",
      author_email="jean.cury@normalesup.org",
      url="https://github.com/gem-pasteur/Integron_Finder/",
      download_url='https://github.com/gem-pasteur/Integron_Finder/archive/v1.5.1.tar.gz',
      license="GPLv3",
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 2',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
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
      data_files=expand_data([('share/integron_finder/data/', ['data']),
                              ('share/integron_finder/doc/html', ['doc/build/html']),
                              ('share/integron_finder/doc/pdf', ['doc/build/latex/IntegronFinder.pdf'])
                             ]),
      cmdclass={'install_scripts': install_scripts},
      distclass=UsageDistribution
      )


