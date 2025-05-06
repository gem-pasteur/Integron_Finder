.. IntegronFinder - Detection of Integron in DNA sequences

.. _install:

************
Installation
************

.. _dependencies:

IntegronFinder dependencies
===========================

IntegronFinder is built with Python >= 3.10, and a few libraries are needed:

- Python >=3.10
- Numpy >=1.26
- Matplotlib >=3.8
- Pandas >=2
- Biopython >=1.82
- colorlog

From version 1.5.1, integron_finder will check and install theses libraries for you.

In addition, IntegronFinder has external dependencies, which have to be
installed prior the use of the program (click to access the corresponding
website).

- `HMMER`_ >=3.1b2
- `INFERNAL`_ >=1.1.2
- `Prodigal`_ >V2.6.2
- `Nextflow`_ (for parallelization)

After installation of these programs, they should be in your ``$PATH`` (*i.e.*
you can type in a terminal ``hmmsearch``, ``cmsearch``, or ``prodigal`` and a
``command not found`` shall not be displayed). If you have them installed
somewhere else, please refer to integron_finder's parameters to give complete path to
IntegronFinder.

.. _installation:

Installation procedure
======================

.. warning::
    When installing a new version (up to 2.0 included) of IntegronFinder,
    do not forget to :ref:`uninstall <uninstallation>` the previous version installed !

.. warning::
    If You upgrading from version prior to 2.0 to 2.0 be careful the python used changed for 3.x.
    The python 2.7 is not supported anymore. So if you installed ``integron_finder`` within a virtualenv
    you need to create a new one based on python3.


From Version 2.0
----------------

System wide installation
""""""""""""""""""""""""
1. Open a terminal and hit (not recommended)::

    sudo pip install integron_finder

.. warning::
    On recent Debian/Ubuntu the --user option is forced. So use of --root option give an unexpected behavior
    and you cannot use --prefix option at all unless you add option --system
    for instance ::

        sudo pip install --system integron_finder

    or ::

        pip install --prefix=/tmp/test_if --system integron_finder


2. To get an updated version (no need to uninstall)::

    sudo pip install -U integron_finder


User wide installation
""""""""""""""""""""""

1. Open a terminal and hit::

    pip install --user integron_finder


Installation in a virtualenv
""""""""""""""""""""""""""""

The virtual environment (`virtualenv`_) is a system to isolate a python program from the system and avoid libraries conflict.
So you can install a different python or libraries version than your system in each virtualenv.
So if you update the system it will not change anything for your program and *vice versa*.
If you want to remove the program just remove the virtual environment.

Create a virtual environment::

    python3 -m venv Integron_Finder

or on some systems::

    virtualenv -p python3 Integron_Finder


activate you virtualenv::

    source Integron_Finder/bin/activate

The name of the virtualenv appear in parenthesis at the beginning of the prompt.
Then install integron_finder::

    python3 -m pip install integron_finder


To run integron finder, you have to activate (once per session) the virtual environment::

    source  Integron_Finder/bin/activate

When you do not need to use integron_finder just deactivate the virtual environment.
In the active terminal just type::

    deactivate

The integron_finder command will disappear from the path.
The name of the virtualenv disappear from the prompt.

Conda Package
"""""""""""""

From 2.0 version, Integron_Finder is available as `conda <https://conda.io/docs/index.html>`_ package.
Integron_finder is in `bioconda <https://bioconda.github.io/>`_ From 2.0 version, Integron_Finder is available as [conda](https://conda.io/docs/index.html) package.
Integron_finder is in [bioconda](https://bioconda.github.io/) channel.
(The advantage with this solution is that it will install prodigal, hmmer, and infernal too.)

1. install conda
2. Set up channels ::

    conda config --add channels conda-forge
    conda config --add channels bioconda

3. install integron_finder ::

    conda install integron_finder

   (The advantage with this solution is that it will install prodigal, hmmer, and infernal too.)


From Version 1.5.1 and after
----------------------------

1. Open a terminal and hit::

    (sudo) pip install integron_finder

2. To get an updated version (no need to uninstall)::

    (sudo) pip install -U integron_finder

For Version 1.5 and before
--------------------------

1. Download the `latest release`_ that can be installed like this (v1.5)
2. Uncompress it
3. In a shell (*e.g.* a terminal), go to the directory and run::

    (sudo) python setup.py install


.. note::
  Super-user privileges (*i.e.*, ``sudo``) are necessary if you want to
  install the program in the general file architecture.

.. note::
  If you do not have the privileges, or if you do not want to install
  IntegronFinder in the Python libraries of your system, you can install
  IntegronFinder in a virtual environment. See `virtualenv`_ or if you're using
  Canopy, see `Canopy CLI`_

.. warning::
    The installer does not work with pure setuptools procedure, it does not work in egg.
    Unless you disable egg by using the ``--root`` option.
    ``python setup.py install --root /prefix/where/to/install/integron_finder``

.. _uninstallation:

Uninstallation procedure
=========================

From Version 1.5.1 and after
----------------------------

To uninstall IntegronFinder, run in the following command::

    (sudo) pip uninstall integron_finder

It will uninstall integron_finder executable

From Version 1.0 to Version 1.5
-------------------------------

Go to the directory from where you installed IntegronFinder (e.g. Integron_Finder-1.5), and run::

    (sudo) python setup.py uninstall

.. _python_install:

How to install Python
=====================

The purpose of this section is to provide some help about installing python
dependencies for IntegronFinder if you never installed any python package.

As IntegronFinder has not been tested on Windows, we assume Unix-based operating system.
For Windows users, the best would be to install a unix virtual machine on your computer.

Usually a python distribution is already installed on your machine.
However, if you don't know how to install libraries, we recommend to re-install it from a distribution which contains pre-compiled libraries.
There are two main distributions (click to access website):

- `Enthought Canopy`_
- `conda-forge`_

Download version 3.x which correspond to your machine,
then make sure that python from these distributions is the default one
(you can possibly choose that in the preference and/or during installation).
Make sure Biopython is installed, otherwise, you will have to install Biopython.
``pip`` or ``conda`` are recommended as a python packages installer.

It works as follow::

    (sudo) pip install Biopython==1.82

To install version 1.82 of Biopython (recommended for IntegronFinder).

.. note::
    If you don't manage to install all the packages, try googling the error, or don't hesitate to ask a question on `stackoverflow`_.

.. _`conda-forge`: https://conda-forge.org/
.. _`Enthought Canopy`: https://store.enthought.com/
.. _`Canopy for Academics`: https://store.enthought.com/#canopy-academic
.. _`stackoverflow`: http://stackoverflow.com/

.. _`HMMER`: http://hmmer.janelia.org/
.. _`INFERNAL`: http://infernal.janelia.org/
.. _`Prodigal`: https://github.com/hyattpd/Prodigal/releases
.. _`Nextflow`: https://www.nextflow.io/

.. _`latest release`: https://github.com/gem-pasteur/Integron_Finder/releases
.. _`virtualenv`: https://docs.python.org/3/library/venv.html
.. _`Canopy CLI`: http://docs.enthought.com/canopy/configure/canopy-cli.html#canopy-cli-venv
