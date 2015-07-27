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
   - Pandas 0.15.1
   - Numpy 1.9.1
   - Biopython 1.63
   - Matplotlib 1.4.3

In addition, IntegronFinder have external dependencies which have to be
installed prior the use of the program.

- `HMMER 3.1b1`_
- `INFERNAL 1.1`_
- `Prodigal V2.6.2`_

After installation of these programs, they should be in your ``$PATH`` (*i.e.*
you can type in a terminal ``hmmsearch``, ``cmsearch``, or ``prodigal`` and a
``command not found`` shall not be displayed). If you have them installed
somewhere else, please refers to the parameters to give complete path to
IntegronFinder.

.. _installation:

Installation procedure
======================

1. Download the `latest release`_
2. Uncompress it
3. In a shell (*e.g.* terminal), go to the directory::

     cd Integron_Finder-x.x/

4. Start installation with::

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
  When installing a new version of IntegronFinder, do not forget to uninstall the previous version installed !

Uninstallation procedure
=========================

To uninstall IntegronFinder, run in the ``Integron_Finder-x.x/`` directory::

    (sudo) python setup.py uninstall


.. _`HMMER 3.1b1`: http://hmmer.janelia.org/
.. _`INFERNAL 1.1`: http://infernal.janelia.org/
.. _`Prodigal V2.6.2`: https://github.com/hyattpd/Prodigal/releases
.. _`latest release`: https://github.com/gem-pasteur/Integron_Finder/releases/latest
.. _`virtualenv`: http://www.virtualenv.org/
.. _`Canopy CLI`: http://docs.enthought.com/canopy/configure/canopy-cli.html#canopy-cli-venv
