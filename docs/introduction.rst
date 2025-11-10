==============================
Introduction
==============================

Welome to PetThermoTools - an open-source Python3 tool for performing phase equilibria calculations using the MELTS or Holland et al. (2018) thermodynamic models.

Important Notice
===========
As of v0.2.48 the package is hosted on PyPI as `petthermotools` (i.e., all lower case). This means that installation and import of PetThermoTools must ALWAYS use the lower case name.

Using MELTS
===========

As a default, this package uses the MELTS algorithms and the alphaMELTS for Python package, developed by Dr Paula Antoshechkina,  to perform the thermodynamic calculations. 
Therefore, it is neccessary for any user to first install the alphaMELTS for Python package. 
This is easily done through our installation notebooks that you can find on this page (https://petthermotools.readthedocs.io/en/latest/Installation/InstallationScript.html).

Using MAGEMin
============

This package also provides the user the ability to perform calculations using the thermodynamic databases of Green et al. (2025) and Weller et al. (2024). This is made possible through 
the recent development and release of MAGEMin (https://github.com/ComputationalThermodynamics; https://doi.org/10.1029/2022GC010427). Our Python code is set up to call the MAGEMin_C functions form julia. 
This opens the possibility to directly compare different thermodynamic models in a way that has not been possible before.

Before the PetThermoTools package can be used to perform MAGEMin calculations, MAGEMin_C (Riel et al. 2022) and MAGEMinCalc must be installed.
As with the alphaMELTS for Python package this can be easily achieved using our installation notebooks (https://petthermotools.readthedocs.io/en/latest/Installation/MAGEMin_PythonCall.html).


Citing PetThermoTools
=======================

At present there is no publication associated with PetThermoTools (although there are plans to release a pre-print in Fall 2023).
We request that any users using PetThermoTools in their research cite the PetThermoTools zenodo repository https://zenodo.org/badge/latestdoi/536711798 and state the 
version that was used.

YouTube Channel
=======================
As of November 2025 we have a YouTube channel where we will post videos covering installation, the background of PetThermoTools, and basic examples. Please check it out (https://www.youtube.com/channel/UCI-2VyZlb2lpScz9vlSn35w)!

Reporting bugs/issues with the code
==============================
No software is free of bugs, particularly when still in the early stages of development! Many of the calculations performed in PetThermoTools have been benchmarked to the results obtained by alternative MELTS software packages (e.g., MELTS for Excel). However, if users spot any bugs, or wish to request features, they should submit an 'issue' on the GitHub page. Alternatively, they can email the development team. In both cases, please upload any files the code is using (e.g. excel, jupyter notebooks) so that I can run your code to see what the issue is!





