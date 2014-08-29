# coding: utf-8

""" oracle, the suppository of wisdom """ 

import os
import re
import subprocess
import sys

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

    if "--skip-si" not in map(str.lower, sys.argv):
        # Check for ifort
        print("Checking for ifort compiler..")
        if os.system("which ifort") > 0:
            raise Exception("No ifort compiler found. Please download the trial "\
                "version from https://software.intel.com/en-us/intel-fortran-compo"\
                "ser-xe-evaluation-options")

        # Install SI
        print("Installing SI..")
        cwd = os.path.join(os.path.dirname(os.path.abspath(
            os.path.expanduser(__file__))), "si/code")
        installer = subprocess.call("make", cwd=cwd, shell=True)
        if installer != 0:
            raise Exception("SI fortran code could not be compiled")

    else:
        print("Skipping SI installation")

setup(name="oracle",
    version=version,
    author="Andrew R. Casey",
    author_email="arc@ast.cam.ac.uk",
    packages=["oracle"],
    url="http://www.github.com/andycasey/oracle/",
    license="MIT",
    description="the suppository of all wisdom",
    long_description=readfile(os.path.join(os.path.dirname(__file__), "README.md")),
    install_requires=readfile(
        os.path.join(os.path.dirname(__file__), "requirements.txt")).split("\n"),
    entry_points={
        "console_scripts": ["oracle = oracle.cli:main"]
    },
    scripts=["si/si_lineform"]
)
