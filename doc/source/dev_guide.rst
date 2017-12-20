.. IntegronFinder - Detection of Integron in DNA sequences

***************
Developer guide
***************

This part is for developers, who want to work on IntegronFinder scripts.


Developer installation
======================

To be able to work on the scripts and test your changes, you can either work with the scripts without installing the package, either install the package in developer mode. If you do not install the package, you always need to be in the root directory of IntegronFinder to launch it using::

    ./integron_finder

If you want to install it in developer mode, run, from the root of IntegronFinder::

    pip install -e .

Then, you will be able to run IntegronFinder using ``integron_finder`` from anywhere.

In any case, you also need to define the ``INTEGRON_HOME`` environment variable, with the path to the root of IntegronFinder directory::

    export INTEGRON_HOME=<path_to_integronfinder_directory>


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


