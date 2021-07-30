.. IntegronFinder - Detection of Integron in DNA sequences


This part is for developers, who want to work on IntegronFinder scripts.


.. _install_dev:

Developer installation
======================

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

It installs the requirements and create a directory in the virtualenv src/integron_finder
and create links in the virtualenv. So ``integron_finder`` is runnable and you can modify the sources and run it again
without to reinstall the project.

.. note::
    `[dev]` allow to install extra dependencies to generate documentation, compute test coverage ...

.. warning::
    Debian/Ubuntu distribution `--user` is the default. So the `--prefix` option does not work
    and the `--root` option has unexpected behavior. Therefore the best solution is to use `--user` or a virtualenv.


Build a new release
===================

#. activate the virtualenv::

    ./Integron_Finder/bin/activate

#. Go to the root of the project::

    cd Integron_Finder/src/Integron_Finder

#. Te build the new release::

    python -m build .

it will create a source *tar.gz* distribution and a *wheel*


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


.. _tests:

Tests
=====

IntegronFinder is provided with unit tests. You can find them in ``tests`` directory.
You can use them to check that your changes did not break the previous features,
and you can update them, and add your own tests for the new features.

Tests are done using unittest.

.. _runtest:

Running tests
-------------

To run the tests -v option is to increase the verbosity of the output::

    python setup.py test

or::

    python tests/run_tests.py -vv

or::

    python tests/run_tests.py -vv tests/test_utils.py

to run specific tests.

If you also want to get code coverage (you need to install coverage)::

    coverage run  --source integron_finder tests/run_tests.py

Add ``-vv`` to get more details on each test passed/failed.
If you want to see the coverage in html output, run (after executing the command above)::

     coverage html

The html coverage report will be generated in ``coverage_html/index.html``.


Adding tests
------------

If you want to create a new test file, adding a file in tests directory, must start with ``test_``.
Then, write your TestCase by inherits from IntegronTest and your tests using unittest framework
(see examples in existing files), and :ref:`run them<runtest>`.


.. _documentation:

Documentation
=============

Documentation is done using ``sphinx``. Source files are located in ``doc/sources``.
To generate the documentation you just have to run the makefile located in *doc* directory. ::

    make html

To generate the documentation in *html* format or ::

    make latexpdf

to generate the documentation in pdf format (for this option you need to have latex installed on your compute)

You can complete them.
