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

**Contributing to ``oracle``**

If you would like to contribute something to ``oracle`` (e.g., code, line lists, test data) then **please do so**. Follow these steps:

1. [Fork the repository](https://github.com/andycasey/oracle/fork). This will create a copy of the repository under your username so
   that you can make any changes you would like. Your forked repository will be at ``github.com/my_username/oracle.git``.

2. Clone your repository into a new directory by using the following:

   ````
   git clone git@github.com:my_username/oracle.git my_forked_oracle_repository
   ````

3. Make the changes you would like to make in the ``my_forked_oracle_repository`` version.

4. From the ``my_forked_oracle_repository`` in a terminal, add the changes:

   ````
   # This command will show you what has changed.
   git status

   # Let's say you created a new file called "extra_code.py" and you removed a file called "old_test_data.txt"
   git add extra_code.py
   git commit -m "Added new code to do super great awesome stuff"
 
   git rm -f old_test_data.txt
   git commit -m "Removed old test data because we don't need it anymore"

   # Then when you're done making changes, push them to your forked repository on GitHub
   git push
   ````

5. Your changes will now be present in your forked repository. To introduce these changes into the main repository (``andycasey/oracle``),
   [create a pull request](https://github.com/andycasey/oracle/compare) describing the changes you have made.

The process might seem long but it's to keep a good version control history, and it means (later on down the track) we can automatically
run tests on new code in a pull request, to make sure that the changes don't break anything. It also cements your contributions to this
work for anyone else to see!


** Reporting bugs and requesting features **
   
If it's a bug, please [open an issue](http://github.com/andycasey/oracle/issues/new) describing what you were trying to do, what you
expected the behaviour to be, and how the code actually behaved. Include any relevant log/terminal output, figures or links to test data.

If it's a feature request, please [open an issue](http://github.com/andycasey/oracle/issues/new) describing the feature you would like.
Be as specific as possible and include any relevant equations/figures.
