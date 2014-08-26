**oracle, the suppository of all knowledge**
============================================

**Getting Started**

1. [This guide](https://help.github.com/articles/set-up-git) will step you through installing git,
   and will show you how to set up ``git`` with your GitHub credentials so that you can access 
   private GitHub repositories like this from the command line.
 
2. Make your life easier by setting up ``git`` with [SSH keys](https://help.github.com/articles/generating-ssh-keys).

3. To use ``SI`` (which ``oracle`` depends on), you will need the Intel Fortran compiler. You can
   download a trial version of it for free from [here](https://software.intel.com/en-us/intel-fortran-composer-xe-evaluation-options).

4. Now you can install ``oracle``:

   ````
   # This will create an 'oracle' directory in your current working directory 
   git clone git@github.com:andycasey/oracle.git 
   cd oracle
   python setup.py install --user
   # OR if you have sudo access:
   sudo python setup.py install 
   ````

**Updating ``oracle``**

If there are changes to ``oracle`` you can update by doing:

````
cd /wherever/oracle/is/kept
# The following line will download the changes from this GitHub repository
git pull
# Now we will install the latest version so that it's available anywhere on the system 
python setup.py install --user
# OR if you have sudo access:
sudo python setup.py install
````

Any problems? [Open an issue](http://github.com/andycasey/oracle/issues/new).
