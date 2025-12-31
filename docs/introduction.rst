==============================
Introduction
==============================

Welome to PetThermoTools - an open-source Python3 tool for performing phase equilibria calculations through alphaMELTS or MAGEMin.

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


Citing PetThermoTools and underlying packages
=======================

PetThermoTools was primarily designed as a research and teaching tool, allowing users to run thermodynamic calculations with a range of different models and approaches while requiring little-to-no prior Python knowledge. This work does not exist in isolation and utilizes some incredible work that has come before, both in developing thermodynamic models and in the construction of minimization routines. We request that if you use PetThermoTools in your research that you cite the corresponding preprint (Gleeson et al. 2025) and **all of the underlying models and software packages** that your study uses in the main text.  

Please use the following guidlines for citing the correct work in the main text of your manuscript:

**For all studies using PetThermoTools:**

Gleeson, M.L.M., Wieser, P.E., Antoshechkina, P. and Riel, N., 2025. PetThermoTools: a fast, ﬂexible, and accessible Python3 package for performing thermodynamic calculations. Earth ArXiv. https://doi.org/10.31223/X5ZX8F

**For all studies using any of the MELTS thermodynamic models please cite**

*The original Adiabat_1ph release (this is the precursor to the current alphaMELTS distribution):*

Smith, P.M. and Asimow, P.D., 2005. Adiabat_1ph: A new public front‐end to the MELTS, pMELTS, and pHMELTS models. Geochemistry, Geophysics, Geosystems, 6(2). https://doi.org/10.1029/2004GC000816 

*The AGU abstract accompanying the original release of MELTS for MATLAB:*

Antoshechkina, P.M. and Ghiorso, M.S., 2018, December. MELTS for MATLAB: A new educational and research tool for computational thermodynamics. In AGU Fall Meeting Abstracts (Vol. 2018, pp. ED44B-23).

*The studies that underpin the MELTS thermodynamic models and algorithms:*

Ghiorso, M.S. and Sack, R.O., 1995. Chemical mass transfer in magmatic processes IV. A revised and internally consistent thermodynamic model for the interpolation and extrapolation of liquid-solid equilibria in magmatic systems at elevated temperatures and pressures. Contributions to Mineralogy and Petrology, 119(2), pp.197-212. https://doi.org/10.1007/BF00307281 

Asimow, P.D. and Ghiorso, M.S., 1998. Algorithmic modifications extending MELTS to calculate subsolidus phase relations. American mineralogist, 83(9-10), pp.1127-1132. https://doi.org/10.2138/am-1998-9-1022

*If you are using the pMELTS thermodynamic model please cite:*

Ghiorso, M.S., Hirschmann, M.M., Reiners, P.W. and Kress III, V.C., 2002. The pMELTS: A revision of MELTS for improved calculation of phase relations and major element partitioning related to partial melting of the mantle to 3 GPa. Geochemistry, Geophysics, Geosystems, 3(5), pp.1-35. https://doi.org/10.1029/2001GC000217 

*If you are using the rhyolite-MELTS v1.0.2 model please cite:*

Gualda, G.A., Ghiorso, M.S., Lemons, R.V. and Carley, T.L., 2012. Rhyolite-MELTS: a modified calibration of MELTS optimized for silica-rich, fluid-bearing magmatic systems. Journal of Petrology, 53(5), pp.875-890. https://doi.org/10.1093/petrology/egr080 

*If you are using the rhyolite-MELTS v1.1.0 or rhyolite-MELTS v1.2.0 models please cite Gualda et al. (2012) and:*

Ghiorso, M.S. and Gualda, G.A., 2015. An H2O–CO2 mixed fluid saturation model compatible with rhyolite-MELTS. Contributions to Mineralogy and Petrology, 169(6), p.53. https://doi.org/10.1007/s00410-015-1141-8

**For all studies using the Green et al. (2025) or Weller et al. (2024) thermodynamic models please cite**

*The recent release of MAGEMin and MAGEMin_C:*

Riel, N., Kaus, B.J., Green, E.C.R. and Berlie, N., 2022. MAGEMin, an efficient Gibbs energy minimizer: application to igneous systems. Geochemistry, Geophysics, Geosystems, 23(7), p.e2022GC010427. https://doi.org/10.1029/2022GC010427 

*If you are using the Green et al. (2025) thermodynamic model please cite:*
Green, E.C., Holland, T.J., Powell, R., Weller, O.M. and Riel, N., 2025. Corrigendum to: Melting of Peridotites through to Granites: a Simple Thermodynamic Model in the System KNCFMASHTOCr, and, a Thermodynamic Model for the Subsolidus Evolution and Melting of Peridotite. Journal of Petrology, 66(1), p.egae079. https://doi.org/10.1093/petrology/egae079

Holland, T.J., Green, E.C. and Powell, R., 2018. Melting of peridotites through to granites: a simple thermodynamic model in the system KNCFMASHTOCr. Journal of Petrology, 59(5), pp.881-900. https://doi.org/10.1093/petrology/egy048

*If you are using the Weller et al. (2024) thermodynamic model please cite:*
Weller, O.M., Holland, T.J., Soderman, C.R., Green, E.C., Powell, R., Beard, C.D. and Riel, N., 2024. New thermodynamic models for anhydrous alkaline-silicate magmatic systems. Journal of Petrology, 65(10), p.egae098. https://doi.org/10.1093/petrology/egae098 

YouTube Channel
=======================
As of November 2025 we have a YouTube channel where we will post videos covering installation, the background of PetThermoTools, and basic examples. Please check it out (https://www.youtube.com/channel/UCI-2VyZlb2lpScz9vlSn35w)!

Reporting bugs/issues with the code
==============================
No software is free of bugs, particularly when still in the early stages of development! Many of the calculations performed in PetThermoTools have been benchmarked to the results obtained by alternative MELTS software packages (e.g., MELTS for Excel). However, if users spot any bugs, or wish to request features, they should submit an 'issue' on the GitHub page. Alternatively, they can email the development team. In both cases, please upload any files the code is using (e.g. excel, jupyter notebooks) so that I can run your code to see what the issue is!





