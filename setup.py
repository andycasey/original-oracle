# coding: utf-8

""" oracle, the suppository of wisdom """ 

import os
import re
import subprocess
import sys
from glob import glob

try:
    from setuptools import setup

except ImportError:
    from distutils.core import setup

major, minor1, minor2, release, serial =  sys.version_info
open_kwargs = {"encoding": "utf-8"} if major >= 3 else {}

def readfile(filename):
    with open(filename, **open_kwargs) as fp:
        contents = fp.read()
    return contents

version_regex = re.compile("__version__ = \"(.*?)\"")
contents = readfile(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "oracle", "__init__.py"))

version = version_regex.findall(contents)[0]

# This is the hackiest thing ever, but it's worked ok for my MOOGSILENT package
# and I don't know how to properly set this up with distutils yet (Higher Fortran
# Mantra required to continue).

if "install" in sys.argv:

    if "skip-si" not in map(str.lower, sys.argv):
        # Check for ifort
        print("Checking for ifort compiler..")
        if subprocess.call("which ifort", shell=True, env=os.environ.copy()) > 0:
            raise Exception("No ifort compiler found. Please download the trial "\
                "version from https://software.intel.com/en-us/intel-fortran-compo"\
                "ser-xe-evaluation-options")

        # Install SI
        cwd = os.path.join(os.path.dirname(os.path.abspath(
            os.path.expanduser(__file__))), "oracle/si/code")
        
        #print("Removing any old instances of SI..")
        #try:
        #    os.remove(os.path.join(cwd, "../si_lineform"))
        #    map(os.remove, glob(os.path.join(cwd, "*.o")))
        #except OSError:
        #    None

        print("Installing SI..")
        installer = subprocess.call("make", cwd=cwd, shell=True,
            env=os.environ.copy())

        # Clean up either way.
        #print("Cleaning up..")
        #cleaner = map(os.remove, glob(os.path.join(cwd, "*.o")))
        
        if installer != 0:
            raise Exception("SI fortran code could not be compiled")

    else:
        print("Skipping SI installation..")
        sys.argv.remove("skip-si")

setup(name="oracle",
    version=version,
    author="Andrew R. Casey",
    author_email="arc@ast.cam.ac.uk",
    packages=["oracle", "oracle.si", "oracle.models"],
    url="http://www.github.com/andycasey/oracle/",
    license="MIT",
    description="the suppository of all wisdom",
    long_description=readfile(os.path.join(os.path.dirname(__file__), "README.md")),
    install_requires=readfile(
        os.path.join(os.path.dirname(__file__), "requirements.txt")).split("\n"),
    entry_points={
        "console_scripts": ["oracle = oracle.cli:main"]
    },
    scripts=["oracle/si/si_lineform"],
    include_package_data=True,
    #package_data={"": [""]}
)
