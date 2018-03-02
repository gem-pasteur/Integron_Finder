.. IntegronFinder - Detection of Integron in DNA sequences


This part is for developers, who want to work on IntegronFinder scripts.


.. _install_dev:

Developer installation
======================

If you are not part of the project, start by forking IntegronFinder repository.
For that, sign in to your account on github, and go to https://github.com/gem-pasteur/Integron_Finder.
Then, click on 'Fork' (under your account icone).
This will create a copy of the repository, but with your username instead of 'gem-pasteur'.

create a virtual environment::

    virtualenv -p python2.7 integron_finder

activate you virtualenv::

    source integron_finder/bin/activate

then install integron_finder in developper mode::

    pip install -e https://github.com/gem-pasteur/Integron_Finder#egg=integron_finder

it install the requirements and create a directory in the virtualenv src/integron_finder
and create links in the virtualenv. So integron finder is runnable and you can modify the sources and run it again
without to reinstall the project.

Send changes to upstream repository
===================================

If you want to integrate your code in the upstream (main) repository, you need to
create a pull request.

1. create a new branch with ``<your branch name>`` a descriptive name
   (e.g. 'adding-xx-feature', 'fixing-typos', etc.), so that others understand what your are working on.
2. work on it
3. test that your work does not break the tests.
   add tests corresponding to your code
4. push your local branch on your integron_finder clone on github ::

   git push --set-upstream origin <your branch name>

4. ask for pull request

    - go to your forked repository on github `https://github.com/<your_login>/Integron_Finder/pulls`
    - Click on 'New pull request'
    - Choose your repository and the branch on which you did your changes in 'head fork' (right-hand side), and choose 'gem-pasteur/Integron_Finder' with the branch on which you want to merge (probably master) in 'base fork' (left-hand side).
    - A green 'Able to merge' text should appear if git is able to automatically merge the 2 branches. In that case, click on 'Create pull request', write your comments on the changes you made, why etc., and save. We will receive the pull request.


.. _tests:

Tests
=====

IntegronFinder is provided with unit tests. You can find them in ``tests`` diretory.
You can use them to check that your changes did not break the previous features,
and you can update them, and add your own tests for the new features.

Tests are done using unittest.

.. _runtest:

Running tests
-------------

To run the tests -v option is to increase the verbosity of the output::

    python tests/run_tests.py -vv

or::

    python tests/run_tests.py -vv tests/test_utils.py

to run specific tests
Or, if you also want to get code coverage::

    coverage run  --source integron_finder tests/run_tests.py

Add ``-vv`` to get more details on each test passed/failed.
If you want to see the coverage in html output, run::

     coverage html

The html coverage report will be generated in ``coverage_html/index.html``.


Adding tests
------------

If you want to create a new test file, adding a file in tests directory, must start with ``test_``.
Then, write your TestCase by inherits from IntegronTest and your tests using unittest framework
(see examples in existing files), and :ref:`run them<runtest>`.

Documentation
=============

Documentation is done using ``sphinx``. Source files are located in ``doc/sources``. You can complete them.

To generate html documentation, go to ``doc`` directory, and run::

    make html


