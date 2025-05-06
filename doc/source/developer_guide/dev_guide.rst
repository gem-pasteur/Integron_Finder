.. IntegronFinder - Detection of Integron in DNA sequences


This part is for developers, who want to work on IntegronFinder scripts.


.. _install_dev:

**********************
Developer installation
**********************

============
Dependencies
============
**Python version >=3.10** is required to run Integron_Finder: https://docs.python.org/3.10/index.html

Integron_Finder has one program dependency:

 - The *HMMER* program, version 3.1 or greater (http://hmmer.org/).
 - The Infernal program version 1.1.2 or greater (http://eddylab.org/infernal/)
 - The Prodigal program version 2.6.2 or greater (https://github.com/hyattpd/Prodigal)

The *hmmsearch, cmsearch and prodigal* programs should be installed (*e.g.*, in the PATH) in order to use integron_finder.
Otherwise, the paths to these executables must be specified via the command-line: see options `--hmmsearch` `--cmsearch` `--prodigal`

Integron_Finder also relies on some Python library dependencies:

 - colorlog
 - pandas
 - matplotlib
 - sphinx
 - sphinx_rtd_theme
 - coverage
 - build
 - ruff
 - pre-commit

These dependencies will be automatically retrieved and installed when using `pip` for installation (see below).


If you are not part of the project, start by forking IntegronFinder repository.
For that, sign in to your account on github, and go to https://github.com/gem-pasteur/Integron_Finder.
Then, click on 'Fork' (under your account icon).
This will create a copy of the repository, but with your username instead of 'gem-pasteur'.

create a virtual environment::

    virtualenv -p python3 Integron_Finder

activate you virtualenv::

    source Integron_Finder/bin/activate

then install integron_finder in developer mode::

    pip install -e "git+https://github.com/gem-pasteur/Integron_Finder#egg=integron_finder[dev]"

or clone your repository manually, then install it ::

    mkdir src
    cd src
    git clone https://github.com/gem-pasteur/Integron_Finder
    cd Integron_Finder
    pip install -e ".[dev]"

install tools to ensure coding style::

    pre-commit install

It installs the requirements and create a directory in the virtualenv src/integron_finder
and create links in the virtualenv. So ``integron_finder`` is runnable and you can modify the sources and run it again
without to reinstall the project.

.. note::
    `[dev]` allow to install extra dependencies to generate documentation, compute test coverage ...

.. note::

    from 2.0.6 version, *Integron_Finder* has adopted `ruff <https://docs.astral.sh/ruff/>`_ as linter
    and *pre-commit* to ensure the coding style.
    please read `CONTRIBUTING.md <https://github.com/gem-pasteur/macsyfinder/blob/master/CONTRIBUTING.md>`_ guide lines.


.. warning::
    Debian/Ubuntu distribution `--user` is the default. So the `--prefix` option does not work
    and the `--root` option has unexpected behavior. Therefore the best solution is to use `--user` or a virtualenv.


.. _tests:

==================================
Integron_Finder testing procedures
==================================

Integron_Finder project use `unittest` framework (included in the standard library) to test the code.

All tests stuff is in `tests` directory.

* The data directory contains data needed by the tests
* in the __init__.py file a IntegronTest class is defined and should be the base of all testcase use in the project
* each test_*.py represent a file containing unit or functional tests

To run all the tests (in the virtualenv)

.. code-block:: shell

    python -m unittest discover

To increase verbosity of output

.. code-block:: shell

    python -m unittest discover -vv

.. code-block:: text

    test_integron_1elem_int (tests.test_add_feature.TestAddFeature.test_integron_1elem_int) ... ok
    test_integron_1elem_prom (tests.test_add_feature.TestAddFeature.test_integron_1elem_prom) ... ok
    test_integron_1elem_prot (tests.test_add_feature.TestAddFeature.test_integron_1elem_prot) ... ok
    test_integron_2int_nelem (tests.test_add_feature.TestAddFeature.test_integron_2int_nelem) ... ok
    test_integron_long_seqname (tests.test_add_feature.TestAddFeature.test_integron_long_seqname) ... ok
    test_integron_nelem (tests.test_add_feature.TestAddFeature.test_integron_nelem) ... ok
    test_config (tests.test_config.TestConfig.test_config) ... ok
    test_default_topology (tests.test_config.TestConfig.test_default_topology) ... ok
    test_func_annot_path (tests.test_config.TestConfig.test_func_annot_path) ... ok
    test_getattr (tests.test_config.TestConfig.test_getattr) ... ok
    test_input_dir (tests.test_config.TestConfig.test_input_dir) ... ok
    test_log_level (tests.test_config.TestConfig.test_log_level) ... ok
    test_model_attc_name (tests.test_config.TestConfig.test_model_attc_name) ... ok
    test_model_attc_path (tests.test_config.TestConfig.test_model_attc_path) ... ok
    test_model_dir (tests.test_config.TestConfig.test_model_dir) ... ok
    test_model_integrase (tests.test_config.TestConfig.test_model_integrase) ... ok
    test_model_len (tests.test_config.TestConfig.test_model_len) ... ok
    test_model_phage_int (tests.test_config.TestConfig.test_model_phage_int) ... ok
    ...

    ----------------------------------------------------------------------
    Ran 246 tests in 400.863s

    OK

The tests must be in python file (`.py`) starting with with `test\_` \
It's possible to specify one or several test files, one module, or one class in a module or a method in a Test class.

Test the `test_pprot_db` module

.. code-block:: shell

    python -m unittest -vv tests.test_prot_db

.. code-block:: text

    test_ProteinDB (tests.test_prot_db.TestCustomDB.test_ProteinDB) ... ok
    test_ProteinDB_bad_parser (tests.test_prot_db.TestCustomDB.test_ProteinDB_bad_parser) ... ok
    test_coding_prot_ids (tests.test_prot_db.TestCustomDB.test_coding_prot_ids) ... ok
    test_get_description (tests.test_prot_db.TestCustomDB.test_get_description) ... ok
    test_get_description_lazy_parser (tests.test_prot_db.TestCustomDB.test_get_description_lazy_parser) ... ok
    test_get_description_stupid_parser (tests.test_prot_db.TestCustomDB.test_get_description_stupid_parser) ... ok
    test_get_description_stupid_parser2 (tests.test_prot_db.TestCustomDB.test_get_description_stupid_parser2) ... ok
    test_getitem (tests.test_prot_db.TestCustomDB.test_getitem) ... ok
    test_iter (tests.test_prot_db.TestCustomDB.test_iter) ... ok
    test_protfile (tests.test_prot_db.TestCustomDB.test_protfile) ... ok
    test_ProteinDB (tests.test_prot_db.TestGemBase.test_ProteinDB) ... ok
    test_codig_prot_ids (tests.test_prot_db.TestGemBase.test_codig_prot_ids) ... ok
    test_find_gembase_file_basename (tests.test_prot_db.TestGemBase.test_find_gembase_file_basename) ... ok
    ...
    test_ProteinDB_no_prodigal (tests.test_prot_db.TestProdigalDB.test_ProteinDB_no_prodigal) ... ok
    test_coding_prot_ids (tests.test_prot_db.TestProdigalDB.test_coding_prot_ids) ... ok
    test_get_description (tests.test_prot_db.TestProdigalDB.test_get_description) ... ok
    test_getitem (tests.test_prot_db.TestProdigalDB.test_getitem) ... ok
    test_iter (tests.test_prot_db.TestProdigalDB.test_iter) ... ok
    test_make_protfile (tests.test_prot_db.TestProdigalDB.test_make_protfile) ... ok
    test_make_protfile_no_dir (tests.test_prot_db.TestProdigalDB.test_make_protfile_no_dir) ... ok
    test_make_protfile_prodigal_failed (tests.test_prot_db.TestProdigalDB.test_make_protfile_prodigal_failed) ... ok
    test_protfile (tests.test_prot_db.TestProdigalDB.test_protfile) ... ok
    test_str (tests.test_prot_db.TestRepliconType.test_str) ... ok

    ----------------------------------------------------------------------
    Ran 41 tests in 5.249s

    OK


Test only the class `TestProdigalDB` (this module contains 5 classes)

.. code-block:: shell

    python -m unittest -vv tests.test_prot_db.TestProdigalDB

.. code-block:: text

    test_ProteinDB (tests.test_prot_db.TestProdigalDB.test_ProteinDB) ... ok
    test_ProteinDB_no_prodigal (tests.test_prot_db.TestProdigalDB.test_ProteinDB_no_prodigal) ... ok
    test_coding_prot_ids (tests.test_prot_db.TestProdigalDB.test_coding_prot_ids) ... ok
    test_get_description (tests.test_prot_db.TestProdigalDB.test_get_description) ... ok
    test_getitem (tests.test_prot_db.TestProdigalDB.test_getitem) ... ok
    test_iter (tests.test_prot_db.TestProdigalDB.test_iter) ... ok
    test_make_protfile (tests.test_prot_db.TestProdigalDB.test_make_protfile) ... ok
    test_make_protfile_no_dir (tests.test_prot_db.TestProdigalDB.test_make_protfile_no_dir) ... ok
    test_make_protfile_prodigal_failed (tests.test_prot_db.TestProdigalDB.test_make_protfile_prodigal_failed) ... ok
    test_protfile (tests.test_prot_db.TestProdigalDB.test_protfile) ... ok

    ----------------------------------------------------------------------
    Ran 10 tests in 0.857s

    OK


Test only the method `test_protfile` from the test Class `TestProdigalDB` in module `test_prot_db`

.. code-block:: shell

    python -m unittest -vv tests.test_prot_db.TestProdigalDB.test_protfile

.. code-block:: text

    test_protfile (tests.test_prot_db.TestProdigalDB.test_protfile) ... ok

    ----------------------------------------------------------------------
    Ran 1 test in 0.112s

    OK


Coverage
========

To compute the tests coverage, we use the `coverage <https://pypi.org/project/coverage/>`_ package.
The package is automatically installed if you have installed `integron_finder` with the `dev` target see :ref:`installation <install_dev>`
The coverage package is setup in the `pyproject.toml` configuration file

To compute the coverage

.. code-block:: shell

    coverage run

.. code-block:: text

    ...
    test_w_chunk (tests.test_split.TestMain.test_w_chunk) ... ok
    test_wo_chunk (tests.test_split.TestMain.test_wo_chunk) ... ok
    test_mute (tests.test_split.TestParseArgs.test_mute) ... ok
    test_parse_chunk (tests.test_split.TestParseArgs.test_parse_chunk) ... ok
    test_parse_outdir (tests.test_split.TestParseArgs.test_parse_outdir) ... ok
    test_parse_replicon (tests.test_split.TestParseArgs.test_parse_replicon) ... ok
    test_quiet (tests.test_split.TestParseArgs.test_quiet) ... ok
    test_verbose (tests.test_split.TestParseArgs.test_verbose) ... ok
    test_split_avoid_overwriting (tests.test_split.TestSplit.test_split_avoid_overwriting) ... ok
    test_split_w_chunk (tests.test_split.TestSplit.test_split_w_chunk) ... ok
    test_split_wo_chunk (tests.test_split.TestSplit.test_split_wo_chunk) ... ok
    test_getitem (tests.test_topology.TestTopology.test_getitem) ... ok
    test_getitem_cmdline_topofile (tests.test_topology.TestTopology.test_getitem_cmdline_topofile) ... ok
    test_getitem_gembase (tests.test_topology.TestTopology.test_getitem_gembase) ... ok
    test_parse (tests.test_topology.TestTopology.test_parse) ... ok
    test_parse_topology (tests.test_topology.TestTopology.test_parse_topology) ... ok
    test_FastaIterator (tests.test_utils.TestUtils.test_FastaIterator) ... ok
    test_FastaIterator_test_topologies (tests.test_utils.TestUtils.test_FastaIterator_test_topologies) ... ok
    test_get_name_from_path (tests.test_utils.TestUtils.test_get_name_from_path) ... ok
    test_log_level (tests.test_utils.TestUtils.test_log_level) ... ok
    test_model_len (tests.test_utils.TestUtils.test_model_len) ... ok
    test_read_multi_prot_fasta (tests.test_utils.TestUtils.test_read_multi_prot_fasta) ... ok

    ----------------------------------------------------------------------
    Ran 246 tests in 400.863s

    OK


Then display a report

.. code-block:: shell

    coverage report


.. code-block:: text

    Name                                  Stmts   Miss Branch BrPart  Cover
    -----------------------------------------------------------------------
    integron_finder/__init__.py              81     12     18      4    84%
    integron_finder/annotation.py            84      0     30      0   100%
    integron_finder/argparse_utils.py        14      1      2      1    88%
    integron_finder/attc.py                 152      0     62      4    98%
    integron_finder/config.py               104      2     48      2    97%
    integron_finder/hmm.py                   89      0     28      1    99%
    integron_finder/infernal.py             147      0     60      2    99%
    integron_finder/integrase.py             33      0     12      2    96%
    integron_finder/integron.py             403     11    146      8    97%
    integron_finder/prot_db.py              374     17    122      7    95%
    integron_finder/results.py               45      0     10      0   100%
    integron_finder/scripts/__init__.py       0      0      0      0   100%
    integron_finder/scripts/finder.py       284     28    116     19    87%
    integron_finder/scripts/merge.py         82      3     30      4    94%
    integron_finder/scripts/split.py         76      5     28      6    89%
    integron_finder/topology.py              39      0     22      1    98%
    integron_finder/utils.py                 96      1     26      2    98%
    -----------------------------------------------------------------------
    TOTAL                                  2103     80    760     63    95%


or generate a html report

.. code-block:: shell

    coverage html

.. code-block:: text

    Wrote HTML report to htmlcov/index.html

The results are in the `htmlcov` directory. With you favourite web browser, open the `index.html` file.
for more options please refer to the `coverage documentation <https://coverage.readthedocs.io/en/latest/>`_ .

.. _documentation:

=============
Documentation
=============

Documentation is done using ``sphinx``. Source files are located in ``doc/sources``.
To generate the documentation you just have to run the makefile located in *doc* directory. ::

    make html

To generate the documentation in *html* format or ::

    make latexpdf

to generate the documentation in pdf format (for this option you need to have latex installed on your compute)

You can complete them.

===================
Build a new release
===================

#. activate the virtualenv::

    ./Integron_Finder/bin/activate

#. Go to the root of the project::

    cd Integron_Finder/src/Integron_Finder

#. Te build the new release::

    python -m build .

it will create a source *tar.gz* distribution and a *wheel*

===================================
Send changes to upstream repository
===================================

If you want to integrate your code in the upstream (main) repository, you need to
create a pull request.

1. Read the `Contibuting guide <https://github.com/gem-pasteur/Integron_Finder/blob/master/CONTRIBUTING.md>`_
2. Create a new branch with ``<your branch name>`` a descriptive name
   (e.g. 'adding-xx-feature', 'fixing-typos', etc.), so that others understand what your are working on.
3. Work on it
4. Test that your work does not break the tests.
   add tests corresponding to your code
5. Push your local branch on your integron_finder clone on github ::

        git push --set-upstream origin <your branch name>

6. ask for pull request

    - Go to your forked repository on github `https://github.com/<your_login>/Integron_Finder/pulls`
    - Click on 'New pull request'
    - Choose your repository and the branch on which you did your changes in 'head fork' (right-hand side),
      and choose 'gem-pasteur/Integron_Finder' with the branch on which you want to merge
      (probably master) in 'base fork' (left-hand side).
    - A green 'Able to merge' text should appear if git is able to automatically merge the 2 branches.
      In that case, click on 'Create pull request', write your comments on the changes you made, why etc,
      and save. We will receive the pull request.
