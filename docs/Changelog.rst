================================================
Change Log
================================================


Version 0.1.25 (August 8th, 2023)
================================
Added both "comp" and "bulk" keyword arguments to multi_path for consistency with prior versions.

Version 0.1.28 (October 11th, 2023)
================================
Added adiabatic decompression melting function

Version 0.1.29 (November 10th, 2023)
================================
Updated function name for mineral barometry code ('find_mineral_cosaturation').
Added schematic showing current available functions to ReadTheDocs.

For full functionality when using MAGEMin calculations, make sure PetThermoTools v0.1.29 is installed:
pyMAGEMINcalc v0.0.7 https://github.com/gleesonm1/pyMAGEMINcalc/archive/refs/tags/v0.0.7.zip
MAGEMinCalc v0.2.22 https://github.com/gleesonm1/MAGEMinCalc/archive/refs/tags/v0.2.22.zip 
MAGEMin_C v1.3.5 https://github.com/ComputationalThermodynamics/MAGEMin_C.jl/archive/refs/tags/v1.3.5.zip 

Version 0.2.1 (November 22nd, 2023)
================================
Repository name updated to PetThermoTools (Petrology Thermodynamics Tools).
Change-log and available functions linked to ReadTheDocs.
pyMAGEMINcalc, MAGEMinCalc, and MAGEMin_C versions are as listed above.


Version 0.2.4 (January 9th, 2024)
================================
Update to crystallization and decompression codes so that the output DataFrames are always the same length.


Version 0.2.5 (January 11th, 2024)
================================
Liquidus code updated to accept single (integer or float) inputs for initial temperature guess.
Path/crysatllization and degassing codes now return cumulate sum of fractionated phases - e.g., Results['Mass']['olivine1_cumsum']