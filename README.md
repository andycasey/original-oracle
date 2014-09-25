**oracle, the suppository of all wisdom**
============================================

Some helpful links:

- [Install guide](https://github.com/andycasey/oracle/wiki/Installation)
- [How to contribute to ``oracle``](https://github.com/andycasey/oracle/wiki/Contributing)
- [Reporting bugs or requesting features](https://github.com/andycasey/oracle/wiki/Reporting-bugs-and-requesting-features)
- [**The current to do list**](https://github.com/andycasey/oracle/wiki/TODO)


At present there are three different solvers (or ways to determine stellar parameters) in ``oracle``. There are three because I wasn't sure which would be best or fastest, so I wrote them all.

**Theremin Solver**

````
oracle solve theremin <configuration.yaml> <spectrum_blue.fits> [, spectrum_red.fits]
````

Performs synthesis on-the-fly and fits abundances to individual atomic lines. Blends are accounted for in each iteration.


**Generative Solver**

````
oracle solve generative <configuration.yaml> <spectrum_blue.fits> [, spectrum_red.fits]
````

This is probably the Right Way™ to do things but needs some work. At present it will be *really* slow because there is no fancy pre-optimisation prior to performing inference.


**Classical Solver**

````
oracle solve classical <configuration.yaml> <spectrum.fits>
````

This performs a classical excitation and ionisation balance using equivalent widths. It's extremely fast but does not take blended lines into account.

To get an idea of some of the flags you can use, do ``oracle solve --help``.





