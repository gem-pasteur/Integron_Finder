#! /usr/bin/env python2
# -*- coding: utf-8 -*-

from __future__ import print_function

import unittest
import sys
import os


def _run(test_files, test_root_path, verbosity=0):
    if not test_files:
        suite = unittest.TestLoader().discover(test_root_path, pattern="test_*.py")
    else:
        test_files = [t for t in test_files if test_root_path in t]
        suite = unittest.TestSuite()
        for test_file in test_files:
            if os.path.exists(test_file):
                if os.path.isfile(test_file):
                    fpath, fname = os.path.split(test_file)
                    suite.addTests(unittest.TestLoader().discover(fpath, pattern=fname))
                elif os.path.isdir(test_file):
                    suite.addTests(unittest.TestLoader().discover(test_file))
            else:
                sys.stderr.write("{0} : no such file or directory\n".format(test_file))

    test_runner = unittest.TextTestRunner(verbosity=verbosity).run(suite)
    return test_runner


def run_unittests(test_files, verbosity=0):
    """
    Execute Unit Tests

    :param test_files: the file names of tests to run.
    of it is empty, discover recursively tests from 'tests/unit' directory.
    a test is python module with the test_*.py pattern
    :type test_files: list of string
    :param verbosity: the verbosity of the output
    :type verbosity: positive int
    :return: True if the test passed successfully, False otherwise.
    :rtype: bool
    """
    test_root_path = os.path.join(os.path.dirname(__file__), 'unit')
    return _run(test_files, test_root_path, verbosity)


def run_functional_tests(test_files, verbosity=0):
    """
    Execute Functional Tests

    :param test_files: the file names of tests to run.
    of it is empty, discover recursively tests from 'tests/unit' directory.
    a test is python module with the test_*.py pattern
    :type test_files: list of string
    :param verbosity: the verbosity of the output
    :type verbosity: positive int
    :return: True if the test passed successfully, False otherwise.
    :rtype: bool
    """
    test_root_path = os.path.join(os.path.dirname(__file__), 'functional')
    return _run(test_files, test_root_path, verbosity)


if __name__ == '__main__':

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("tests",
                        nargs='*',
                        default=False,
                        help="name of test to execute")

    parser.add_argument("--unit",
                        dest='unit',
                        action='store_true',
                        default=False,
                        help="execute unit tests")

    parser.add_argument("--functional",
                        dest='functional',
                        action='store_true',
                        default=False,
                        help="execute functional tests")

    parser.add_argument("-v", "--verbose",
                        dest="verbosity",
                        action="count",
                        help="set the verbosity level of output",
                        default=0
                        )

    args = parser.parse_args()
    if not any((args.unit, args.functional)):
        args.unit, args.functional = True, True

    result_all_tests = []

    old_path = sys.path
    if 'INTEGRON_HOME' in os.environ and os.environ['INTEGRON_HOME']:
        INTEGRON_HOME = os.environ['INTEGRON_HOME']
        if INTEGRON_HOME not in sys.path:
            sys.path.insert(0, INTEGRON_HOME)
    else:
        # the test directory is automatically added by python
        # in sys.path as first element
        pass

    if args.unit:
        print("\n", "#" * 70, sep='')
        print("Test Runner: Unit tests")
        print("#" * 70)

        ############## WORKAROUND ##################
        # integron_finder is a script not a lib
        # to test each function in integron_finder,
        # we need to import integron_finder
        # so we need to transform the script in lib
        # and ad it to the path somewhere a user can write
        # the fisrt element of path is integron_finder.tests

        if not 'INTEGRON_HOME' in os.environ:
            INTEGRON_HOME = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
            sys.path.append(INTEGRON_HOME)
            from tests import which
            integron_finder_script = which('integron_finder')
            if integron_finder_script is None:
                raise RuntimeError('Cannot find integron_finder, do you set INTEGRON_HOME')
        else:
            INTEGRON_HOME = os.environ['INTEGRON_HOME']
            integron_finder_script = os.path.join(INTEGRON_HOME, 'integron_finder')

        fake_lib = os.path.join(INTEGRON_HOME, 'tests', 'fake_lib')
        integron_finder_lib = os.path.join(fake_lib, 'integron_finder.py')

        if not os.path.exists(fake_lib):
            os.mkdir(fake_lib)
        if os.path.exists(integron_finder_lib):
            os.unlink(integron_finder_lib)

        if not os.path.exists(integron_finder_lib):
            print(integron_finder_script)
            print(integron_finder_lib)
            os.symlink(integron_finder_script, integron_finder_lib)

        sys.path.insert(0, fake_lib)
        sys.path.insert(0, INTEGRON_HOME)
        ############ END WORKAROUND ################

        test_runner = run_unittests(args.tests, verbosity=args.verbosity)
        unit_results = test_runner.wasSuccessful()
        result_all_tests.append(unit_results)

    if args.functional:
        print("\n", "#" * 70, sep='')
        print("Test Runner: Functional tests")
        print("#" * 70)

        test_runner = run_functional_tests(args.tests, verbosity=args.verbosity)
        functional_results = test_runner.wasSuccessful()
        result_all_tests.append(functional_results)

    sys.path = old_path
    if all(result_all_tests):
        sys.exit(0)
    else:
        sys.exit(1)
