# ReadTheDocs

For full documentation, installation instructions, and example calculations using PetThermoTools please visit our readthedocs page.

https://PetThermoTools.readthedocs.io/en/latest/

**Note on Reproducibility:** To comply with the open-source spirit of PetThermoTools and underlying tools, any publication should include a clear explanation of your model settings (Pressure, Temperature, and $fO_2$ buffers, etc.) to ensure results can be independently verified. This is in addition to the citation information listed below.

# Citing PetThermoTools and underlying packages
PetThermoTools was primarily designed as a research and teaching tool, allowing users to run thermodynamic calculations with a range of different models and approaches while requiring little-to-no prior Python knowledge. This work does not exist in isolation and utilizes some incredible work that has come before, both in developing thermodynamic models and in the construction of minimization routines. We request that if you use PetThermoTools in your research that you cite the corresponding preprint (Gleeson et al. 2025) and **all of the underlying models and software packages** that your study uses in the main text.  

Please use the following guidlines for citing the correct work in the main text of your manuscript:

#### For all studies using PetThermoTools:
Gleeson, M.L.M., Wieser, P.E., Antoshechkina, P. and Riel, N., 2025. PetThermoTools: a fast, ﬂexible, and accessible Python3 package for performing thermodynamic calculations. Earth ArXiv. https://doi.org/10.31223/X5ZX8F

#### For all studies using any of the MELTS thermodynamic models please cite
*The original Adiabat_1ph release (this is the precursor to the current alphaMELTS distribution):*
Smith, P.M. and Asimow, P.D., 2005. Adiabat_1ph: A new public front‐end to the MELTS, pMELTS, and pHMELTS models. Geochemistry, Geophysics, Geosystems, 6(2). https://doi.org/10.1029/2004GC000816 

*The current distribution of alphaMELTS:*
Antoshechkina, P. and P. Asimow (2025). magma-source/alphaMELTS: Builds for all OS and chip combinations. Version v2.3.2-beta.0. https://doi.org/10.5281/zenodo.17861695.

*The studies that underpin the MELTS thermodynamic models and algorithms:*
Ghiorso, M.S. and Sack, R.O., 1995. Chemical mass transfer in magmatic processes IV. A revised and internally consistent thermodynamic model for the interpolation and extrapolation of liquid-solid equilibria in magmatic systems at elevated temperatures and pressures. Contributions to Mineralogy and Petrology, 119(2), pp.197-212. https://doi.org/10.1007/BF00307281 

Asimow, P.D. and Ghiorso, M.S., 1998. Algorithmic modifications extending MELTS to calculate subsolidus phase relations. American mineralogist, 83(9-10), pp.1127-1132. https://doi.org/10.2138/am-1998-9-1022

**If you are using the pMELTS thermodynamic model please cite:**
Ghiorso, M.S., Hirschmann, M.M., Reiners, P.W. and Kress III, V.C., 2002. The pMELTS: A revision of MELTS for improved calculation of phase relations and major element partitioning related to partial melting of the mantle to 3 GPa. Geochemistry, Geophysics, Geosystems, 3(5), pp.1-35. https://doi.org/10.1029/2001GC000217 

**If you are using the rhyolite-MELTS v1.0.2 model please cite:**
Gualda, G.A., Ghiorso, M.S., Lemons, R.V. and Carley, T.L., 2012. Rhyolite-MELTS: a modified calibration of MELTS optimized for silica-rich, fluid-bearing magmatic systems. Journal of Petrology, 53(5), pp.875-890. https://doi.org/10.1093/petrology/egr080 

**If you are using the rhyolite-MELTS v1.1.0 or rhyolite-MELTS v1.2.0 models please cite Gualda et al. (2012) and:**
Ghiorso, M.S. and Gualda, G.A., 2015. An H2O–CO2 mixed fluid saturation model compatible with rhyolite-MELTS. Contributions to Mineralogy and Petrology, 169(6), p.53. https://doi.org/10.1007/s00410-015-1141-8

#### For all studies using the Green et al. (2025) or Weller et al. (2024) thermodynamic models please cite
*The recent release of MAGEMin and MAGEMin_C:*
Riel, N., Kaus, B.J., Green, E.C.R. and Berlie, N., 2022. MAGEMin, an efficient Gibbs energy minimizer: application to igneous systems. Geochemistry, Geophysics, Geosystems, 23(7), p.e2022GC010427. https://doi.org/10.1029/2022GC010427 

**If you are using the Green et al. (2025) thermodynamic model please cite:**
Green, E.C., Holland, T.J., Powell, R., Weller, O.M. and Riel, N., 2025. Corrigendum to: Melting of Peridotites through to Granites: a Simple Thermodynamic Model in the System KNCFMASHTOCr, and, a Thermodynamic Model for the Subsolidus Evolution and Melting of Peridotite. Journal of Petrology, 66(1), p.egae079. https://doi.org/10.1093/petrology/egae079

Holland, T.J., Green, E.C. and Powell, R., 2018. Melting of peridotites through to granites: a simple thermodynamic model in the system KNCFMASHTOCr. Journal of Petrology, 59(5), pp.881-900. https://doi.org/10.1093/petrology/egy048

**If you are using the Weller et al. (2024) thermodynamic model please cite:**
Weller, O.M., Holland, T.J., Soderman, C.R., Green, E.C., Powell, R., Beard, C.D. and Riel, N., 2024. New thermodynamic models for anhydrous alkaline-silicate magmatic systems. Journal of Petrology, 65(10), p.egae098. https://doi.org/10.1093/petrology/egae098 


