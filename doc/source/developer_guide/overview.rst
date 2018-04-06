.. IntegronFinder - Detection of Integron in DNA sequences

.. _overview:

*********************
Architecture Overview
*********************


Project files and directories
=============================

files
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

directories
-----------

integron_finder
    The core of the projects contains integron_finder library
    The **scripts/finder** contain the main entry point.

tests
    Contains all needed for tests, the tests themselves, are a the top level and the
    name must start by ``test_``.
    The data directory contains all data needed to perform the tests.

doc
    Contains the documentation write in sphinx. The **source** directory contains the rst files,
    whereas the **build** directory contains the generated documentation

Singularity
    Contains the definition file for singularity container.

data
    TODO

dist
    This directory is generate when a distribution is created (python setup.py sdist).
