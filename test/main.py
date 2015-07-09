#! /usr/bin/env python
# -*- coding: utf-8 -*-


import os
import sys
import unittest


def run(lib, tests, verbosity = 0):
    """
    run the unit tests and print the results on stderr
    
    :param lib: the path where is the integron_finder
    :type lib: string
    :param tests: the name of test to run. if tests is empty list, discover recursively tests form this directory.
                  a test is python module with the test_*.py pattern
    :type tests: list of string
    :param verbosity: the verbosity of the output
    :type verbosity: int
    """
    if lib not in sys.path:
        sys.path.insert(0, lib)
    if not tests:
        suite = unittest.TestLoader().discover(os.path.dirname(__file__), pattern="test_*.py" ) 
    else:
        suite = unittest.TestSuite()
        for test in tests: 
            if os.path.exists(test):
                if os.path.isfile(test):
                    fpath, fname =  os.path.split( test )
                    suite.addTests(unittest.TestLoader().discover(fpath , pattern = fname )) 
                elif os.path.isdir(test):  
                    suite.addTests(unittest.TestLoader().discover(test)) 
            else:
                sys.stderr.write(  "{0} : no such file or directory\n".format(test)) 

    res = unittest.TextTestRunner(verbosity = verbosity).run(suite)
    return res



if __name__ == '__main__':

    if 'INTEGRON_HOME' in os.environ:
        INTEGRON_HOME = os.environ['INTEGRON_HOME']
    else:
        INTEGRON_HOME = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("tests",
                        nargs = '*',
                        default = [],
                        help = "name of test to execute")

    parser.add_argument("-v", "--verbose" , 
                        dest= "verbosity" , 
                        action="count" , 
                        help= "set the verbosity level of output",
                        default = 0
                        )

    args = parser.parse_args()
    res = run(INTEGRON_HOME, args.tests, args.verbosity)
    