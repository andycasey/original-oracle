# coding: utf-8

""" oracle, the suppository of knowledge """ 

import os
import re
import sys

from signal import alarm, signal, SIGALRM, SIGKILL
from subprocess import PIPE, Popen

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
    # Check for ifort
    print("Checking for ifort compiler..")
    if os.system("which ifort") > 0:
        raise Exception("No ifort compiler found. Please download the trial "\
            "version from https://software.intel.com/en-us/intel-fortran-compo"\
            "ser-xe-evaluation-options")

    timeout = 240
    class Alarm(Exception):
        pass

    def alarm_handler(signum, frame):
        raise Alarm

    cwd = os.path.join(os.path.dirname(os.path.abspath(
        os.path.expanduser(__file__))), "si/code")
    env = {
        "PATH": os.getenv("PATH"),
        "EXEDIR": os.path.join(os.path.dirname(os.path.abspath(
            os.path.expanduser(__file__))), "si"),
    }

    print("Compiling SI..")
    p = Popen(["make"],
        shell=False, bufsize=2056, stdin=PIPE, stdout=PIPE, stderr=PIPE,
        cwd=cwd, env=env, close_fds=True)

    if timeout != -1:
        signal(SIGALRM, alarm_handler)
        alarm(timeout)

    try:
        stdout, stderr = p.communicate()
        print("SI Installer output (code {0}):\n{1}".format(p.returncode, stdout))
        if timeout != -1:
            alarm(0)

    except Alarm:
        # Process might have died before getting to this line so wrap it to
        # avoid "OSError: no such process"
        try:
            os.kill(p.pid, SIGKILL)
        except OSError:
            pass
        raise Exception("The SI install process was killed due to timeout.")

    if p.returncode != 0:
        raise Exception("Could not install SI Fortran code (error {0}):\n{1}\n{2}".format(
            p.returncode, stdout, stderr))

setup(name="oracle",
    version=version,
    author="Andrew R. Casey",
    author_email="arc@ast.cam.ac.uk",
    packages=["oracle"],
    url="http://www.github.com/andycasey/oracle/",
    license="MIT",
    description="the suppository of all knowledge",
    long_description=readfile(os.path.join(os.path.dirname(__file__), "README.md")),
    install_requires=readfile(
        os.path.join(os.path.dirname(__file__), "requirements.txt")).split("\n"),
    entry_points={
        "console_scripts": ["oracle = oracle.cli:main"]
    },
    scripts=["fortran/si_lineform"]
)
