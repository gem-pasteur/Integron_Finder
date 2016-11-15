.. IntegronFinder - Detection of Integron in DNA sequences

.. _install:

************
Installation
************

.. _dependencies:

IntegronFinder dependencies
===========================

IntegronFinder is built with Python 2.7, and a few libraries are needed:

- Python 2.7
   - Pandas (>=0.18.0)
   - Numpy (>=1.9.1)
   - Biopython (>=1.65)
   - Matplotlib (>=1.4.3)

From version 1.5.1, integron_finder will check and install theses libraries for you.

In addition, IntegronFinder has external dependencies, which have to be
installed prior the use of the program (click to access the corresponding
website).

- `HMMER 3.1b1`_
- `INFERNAL 1.1`_
- `Prodigal V2.6.2`_

After installation of these programs, they should be in your ``$PATH`` (*i.e.*
you can type in a terminal ``hmmsearch``, ``cmsearch``, or ``prodigal`` and a
``command not found`` shall not be displayed). If you have them installed
somewhere else, please refer to integron_finder's parameters to give complete path to
IntegronFinder.

.. _installation:

Installation procedure
======================

.. warning::
    When installing a new version (up to 1.5.1 included) of IntegronFinder, do not forget to :ref:`uninstall <uninstallation>` the previous version installed !

From Version 1.5.1 and after
----------------------------

1. Open a terminal and hit::

    (sudo) pip install integron_finder

2. To get a updated version (no need to uninstall)::

    (sudo) pip install -U integron_finder

For Version 1.5 and before
--------------------------

1. Download the `latest release`_
2. Uncompress it
3. In a shell (*e.g.* a terminal), go to the directory and run::

    (sudo) python setup.py install


.. note::
  Super-user privileges (*i.e.*, ``sudo``) are necesserary if you want to
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

As IntegronFinder has not been test on Windows, we assume Unix-based operating system. For Windows users, the best would be to install a unix virtual machine on your computer.

Usually a python distribution is already installed on your machine. However, if you don't know how to install libraries, we recommend to re-install it from a distribution which contains pre-compiled libraries. There are two main distributions (click to access website):

- `Enthought Canopy`_
- `Anaconda`_

Download version 2.7 which correspond to your machine, then make sure that python from these distributions is the default one (you can possibly choose that in the preference and/or during installation).
They both come with all the needed packages but Biopython. If you have a **student email adress** from a university-delivering degree, you can request an academic licence to *Enthough Canopy* (see `Canopy for Academics`_) which will allow you to download additional packages including Biopython.

Otherwise, you will have to install Biopython manually. ``pip`` is recommended as a python packages installer. It works as follow::

    (sudo) pip install Biopython==1.65

To install version 1.65 of Biopython (recommended for IntegronFinder).

.. note::
    If you don't manage to install all the packages, try googling the error, or don't hesisate to ask a question on `stackoverflow`_.

.. _`Enthought Canopy`: https://store.enthought.com/
.. _`Anaconda`: https://www.continuum.io/downloads
.. _`Canopy for Academics`: https://store.enthought.com/#canopy-academic
.. _`stackoverflow`: http://stackoverflow.com/

.. _`HMMER 3.1b1`: http://hmmer.janelia.org/
.. _`INFERNAL 1.1`: http://infernal.janelia.org/
.. _`Prodigal V2.6.2`: https://github.com/hyattpd/Prodigal/releases
.. _`latest release`: https://github.com/gem-pasteur/Integron_Finder/releases/latest
.. _`virtualenv`: http://www.virtualenv.org/
.. _`Canopy CLI`: http://docs.enthought.com/canopy/configure/canopy-cli.html#canopy-cli-venv
