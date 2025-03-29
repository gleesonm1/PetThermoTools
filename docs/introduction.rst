==============================
Introduction
==============================

Welome to PetThermoTools - an open-source Python3 tool for performing phase equilibria calculations using the MELTS or Holland et al. (2018) thermodynamic models.

Using MELTS
===========

As a default, this package uses the MELTS algorithms and the alphaMELTS for Python package, developed by Dr Paula Antoshechkina, 
to perform the thermodynamic calculations. Therefore, it is neccessary for any user to first download the alphaMELTS for Python package. 
Please see the installation page for full details.

Using MAGEMin
============

This package also provides the user to perform certain calculations using the thermodynamic database of Holland et al. (2018). This is made possible through 
the recent development and release of MAGEMin (https://github.com/ComputationalThermodynamics; https://doi.org/10.1029/2022GC010427), with 
Julia functions called from Python used to run the calculations. This opens the possibility to directly compare different thermodynamic models in a way that has not been possible before.

Before the PetThermoTools package can be used to perform MAGEMin calculations, the Julia code used to compile and run MAGEMin, which is hosted in a separate 
repository, must be imported and added to Julia. Please see the installation page for full details.


Citing PetThermoTools
=======================

At present there is no publication associated with PetThermoTools (although there are plans to release a pre-print in Fall 2023).
We request that any users using PetThermoTools in their research cite the PetThermoTools zenodo repository https://zenodo.org/badge/latestdoi/536711798 and state the 
version that was used.

Reporting bugs/issues with the code
==============================
No software is free of bugs, particularly when still in the early stages of development! Many of the calculations performed in PetThermoTools have been benchmarked to the results obtained by alternative MELTS software packages (e.g., MELTS for Excel). However, if users spot any bugs, or wish to request features, they should submit an 'issue' on the GitHub page. Alternatively, they can email the development team. In both cases, please upload any files the code is using (e.g. excel, jupyter notebooks) so that I can run your code to see what the issue is!





