#!/usr/bin/env python
"""
Unit tests

Refer to https://trac.sdss3.org/wiki/Ops/Validation for details.
"""
import os
import sys
import glob
import subprocess

global TestPath
TestPath = None

def setup():
    """Add python directory of this package to PYTHONPATH in os.environ
    and set TestPath
    """
    global TestPath
    TestPath = os.path.dirname(__file__)
    pythonPath = os.path.join(os.path.dirname(os.path.abspath(TestPath)), "python")
    if "PYTHONPATH" in os.environ:
        pythonPath += ":" + os.environ["PYTHONPATH"]
    os.environ["PYTHONPATH"] = pythonPath

def runTests(testList):
    """Run the specified tests.
    
    If this test runner file is included in the list it will be skipped.
    """
    for test in testList:
        if test.endswith(__file__):
            continue
        print "\n*** Running %s ***\n" % (test,)
        subprocess.call(["python", test], env=os.environ)

if __name__ == "__main__":
    setup()
    if len(sys.argv) == 1:
        testList = glob.glob(os.path.join(TestPath, "*.py"))
    else:
        testList = sys.argv[1:]
    runTests(testList)
    
