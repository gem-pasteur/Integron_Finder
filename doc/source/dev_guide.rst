.. IntegronFinder - Detection of Integron in DNA sequences

***************
Developer guide
***************

This part is for developers, who want to work on IntegronFinder scripts.


Getting a local working copy of the repository
==============================================

If you are not part of the project, start by forking IntegronFinder repository. For that, sign in to your account on github, and go to https://github.com/gem-pasteur/Integron_Finder. Then, click on 'Fork' (under your account icone). This will create a copy of the repository, but with your username instead of 'gem-pasteur'.
Then, clone it to your computer to have a local working copy, by opening a terminal window, and typing::

    git clone https://github.com/<your_login>/Integron_Finder.git

or, if you linked your ssh key::

    git clone git@github.com:<your_login>/Integron_Finder.git

Go to the ``Integron_Finder`` repository created: ``cd Integron_Finder`` (referred here after as 'the root directory of IntegronFinder'). Then, it is better if you work on your own branch. For that, type::

    git checkout -b <your branch name>

with ``<your branch name>`` a descriptive name (e.g. 'adding-xx-feature', 'fixing-typos', etc.), so that others understand what your are working on.

You can then :ref:`install <install_dev>` IntegronFinder as explained above, and do your modifications, with their associated :ref:`tests <tests>`. Don't forget to commit your changes at each step.

Send changes to upstream repository
===================================

Once you are done, you can push your changes to your forked repository::

    git push --set-upstream origin <your branch name>

Then, create a pull request, to add your changes to the official repository. For that:

- go to your forked repository on github `https://github.com/<your_login>/Integron_Finder/pulls`
- Click on 'New pull request'
- Choose your repository and the branch on which you did your changes in 'head fork' (right-hand side), and choose 'gem-pasteur/Integron_Finder' with the branch on which you want to merge (probably master) in 'base fork' (left-hand side).
- A green 'Able to merge' text should appear if git is able to automatically merge the 2 branches. In that case, click on 'Create pull request', write your comments on the changes you made, why etc., and save. We will receive the pull request.

.. warning:: Before submitting your pull request, please be sure that you provide the unit tests corresponding to the new features you added.



.. _install_dev:

Developer installation
======================

To be able to work on the scripts and test your changes, you can either work with the scripts without installing the package, either install the package in developer mode. If you do not install the package, you always need to be in the root directory of IntegronFinder to launch it using::

    ./integron_finder

If you want to install it in developer mode, run, from the root of IntegronFinder::

    pip install -e .

Then, you will be able to run IntegronFinder using ``integron_finder`` from anywhere.

In any case, you also need to define the ``INTEGRON_HOME`` environment variable, with the path to the root of IntegronFinder directory::

    export INTEGRON_HOME=<path_to_integronfinder_directory>


.. _tests:

Tests
=====

IntegronFinder is provided with unit (resp. functional) tests. You can find them in ``tests/unit`` (resp. ``tests/functional``) repository. You can use them to check that your changes did not break the previous features, and you can update them, and add your own tests for the new features.

Tests are done using unittest.

.. _runtest:

Running tests
-------------

To run the tests (use option corresponding to unit or functional test according to what you want)::

    python tests/run_tests.py [--unit] [--functional]

Or, if you also want to get code coverage::

    coverage run tests/run_tests.py [--unit] [--functional]

Add ``-vvv`` to get more details on each test passed/failed. The html coverage report will be generated in ``coverage_html/index.html``.


Adding tests
------------

If you want to create a new test file, it must be saved in ``tests/unit`` or ``tests/functional`` according to the type of test, and the test file must start with ``test_``. Then, write your tests using unittest framework (examples in existing files), and :ref:`run them<runtest>`.

Documentation
=============

Documentation is done using ``sphinx``. Source files are located in ``doc/sources``. You can complete them.

To generate html documentation, go to ``doc`` directory, and run::

    make html


