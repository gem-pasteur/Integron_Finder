.. IntegronFinder - Detection of Integron in DNA sequences

.. _overview:

*********************
Architecture Overview
*********************


Project files and directories
=============================


Files
-----

COPYING
    The integron_finder licensing.

COPYRIGHT
    The integron finder copy rights holders.

MANIFEST.in
    What must be or should not included in the distribution.

README.md
    The file to red in first.

requirements.txt
    The requirements need to use integron_finder.

requirements_dev.txt
    The extra requirements to develop on integron_finder.

setup.cfg
    The setup.py configuration file.

setup.py
    The file to define how to build/install/release/test/... integron finder.


Directories
-----------

integron_finder
    The core of the projects contains integron_finder library
    The **scripts/finder** contain the main entry point.

tests
    Contains all needed for tests, the tests themselves, are a the top level and the
    name must start by ``test_``.
    The data directory contains all data needed to perform the tests. (see :ref:`tests` for further details)

doc
    Contains the documentation write in sphinx. The **source** directory contains the .rst files,
    whereas the **build** directory contains the generated documentation.
    To know how to contribute or generate documentation see :ref:`documentation`

Singularity
    Contains the definition file for singularity container.

data
    TODO

dist
    This directory is generated when a distribution is created (``python setup.py sdist``).


Technical overview
==================

The main entry point is in integron_finder/scripts/finder.py
there are 3 functions

:func:`intgeron_finder.scripts.main` which is the real main entry point

``main`` call :func:`scripts/finder.parse_args` which parse the commandline and
generate a :class:`config.Config` object.
and do a loop over replicon and run :func:`intgeron_finder.scripts/find_integron_in_one_replicon`

all results are store in a directory named ``Results_Integron_Finder_<replicon_file_name>`` this directory is created by
:func:`intgeron_finder.scripts/find_integron_in_one_replicon` store results in this directory
or in a subdirectory call tmp_<replicon_id> these subdirectories will be keep only if ``--keep-tmp`` option
is set, otherwise they are removed at the end of the :func:`intgeron_finder.scripts/find_integron_in_one_replicon`

when all replicons are computed the ``main`` function call :func:`integron_finder.utils.merge_results` to gather
all results files ``<replicons_id>.integtrons`` and generate a unique file with these information.

to have details on ``find_integron_in_one_replicon`` works see :ref:`introduction`


