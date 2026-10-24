================================================
Change Log
================================================

Version 1.2.3 (Oct 2026)
====================================
multi_path, equilibrate_multi, findLiq_multi and saturation_pressure now use a more robust calculation engine by default. The previous implementation could reuse a MELTS engine's internal state between unrelated rows of a batch, which was found to occasionally produce silently incorrect results (confirmed in roughly 1 in 8 runs at the scale tested for multi_path); every row now always starts from a freshly initialised engine unless you explicitly ask for reuse.
As part of this change, a stuck or crashed calculation within a batch is now handled individually, instead of the whole batch/round waiting on or being discarded by it, and stray worker processes left behind after a crashed or interrupted calculation are now properly cleaned up rather than lingering in the background.
multi_path now reports the timeout being used at the start of a calculation, and gives a clearer breakdown at the end (completed / returned partial results after a crash / returned partial results after timing out / returned no results) instead of a single success/fail count that could under-report calculations that actually returned usable partial results.
Fixed a bug where multi_path, equilibrate_multi, findLiq_multi and saturation_pressure could silently produce incorrect or incomplete results when a pandas Series (e.g. a DataFrame column, such as T_C=data['Temperature']) was passed directly for T_C/P_bar, rather than a plain number, list, or numpy array.
Faster result assembly for large batch calculations.

Version 1.2.0 (Oct 2026)
====================================
New quiet_engine keyword (default True) on multi_path, the crystallisation/decompression wrapper functions, equilibrate_multi, findCO2_multi, findLiq_multi, supCalc, mineral_cosaturation, saturation_pressure, phaseDiagram_calc and AdiabaticDecompressionMelting: suppresses the MELTS engine's console output for a sizeable speed-up in batch calculations. Set quiet_engine=False to see the raw engine output for debugging.
Faster package import: the nGibbs/torch dependency is now only loaded the first time an nGibbs model is actually used.
Repeated import banner messages (e.g. "alphaMELTS for Python files successfully located") no longer print once per worker process during a batch calculation.
Fixed a crash (IndexError) that could occur in later retry rounds of large phaseDiagram_calc calculations.
Fixed the Fe2O3 molar mass used throughout the package (159.59 -> 159.69 g/mol); results shift slightly, mostly by less than 1e-3 relative.
Progress bars no longer require ipywidgets to be installed.
Fixed a Julia package precompilation failure (WeakRefStrings/Parsers version conflict) that could occur in ptt.install_MAGEMinCalc() / ptt.update_MAGEMinCalc(); installation failures now raise immediately rather than only showing up later when the environment is first used.
findSatPressure and findSatPressure_multi now raise a DeprecationWarning pointing to saturation_pressure().

Version 1.1.9 (Sept 2026)
====================================
Fixed bug that was appearing for fluid saturated models.

Version 1.1.6 (Aug 2026)
====================================
Updated installation script to automatically download alphaMELTS v2.3.2 files.

Version 1.1.4 (Jun 2026)
====================================
Improved version control for julia installation

Version 1.1.3 (Jun 2026)
====================================
Change to Python-Julia conection to allow better integration when performing multiple MAGEMin calculations.

Version 1.1.2 (Jun 2026)
====================================
Pre-release of nGibbs included in package for Volcali workshop

Version 1.0.0 (Jun 2026)
====================================
Full release to accompany publication of manuscript in Volcanica
Includes fix to potential temperature search when using pMELTS (can set Tp as input to melting calcs)

Version 0.4.0 (March 2026)
====================================
Major update to MAGEMinCalc. Please run ptt.update_MAGEMinCalc().
Change to how calculations in MAGEMin are called. 
Ability to specify a starting temperature (rather than a potential temperatue for melting calculations).

Version 0.3.22 (March 2026)
====================================
Update to MAGEMinCalc. Please run ptt.update_MAGEMinCalc()
Fix to minor bug with equilibrium crystallization calcs in MAGEMinCalc.

Version 0.3.22 (March 2026)
====================================
Major update to the way path calculations are computed, allowing the results of partial calculations to be returned. Influences phase diagram and barometry calculations as well. 
Should also lead to a small speed improvement due to a move away from Pandas
Code compatible with Pandas v3.

Version 0.3.21 (February 2026)
====================================
Update to harker plotting function to fix bug with individual models

Version 0.3.20 (January 2026)
====================================
Automatic detection of parallelization based on number of performance cores (rather than total logical processors).
Added examples for saving results to ReadTheDocs.

Version 0.3.16 (January 2026)
====================================
Update to import messages and identification of potential memory issues when using julia

Version 0.3.12 (January 2026)
====================================
Relabelling of T_Liq to T_Liq_C in saturation pressure outputs

Version 0.3.11 (January 2026)
====================================
Fix the no solution issue for solubility calculations

Version 0.3.10 (January 2026)
====================================
Updated installation routine for julia packages

Version 0.3.9 (January 2026)
====================================
Replace zeros with nans in saturation pressure calculations

Version 0.3.8 (January 2026)
====================================
Implement of new installation routine for MAGEMinCalc via Victor.

Version 0.3.7 (January 2026)
====================================
Small update to PyPI distribution to fix error in distribution of v0.3.6

Version 0.3.6 (January 2026)
====================================
Added new plotting options for phase abundance (ptt.phase_mass_comparison())

Version 0.3.4 (January 2026)
====================================
Fix minor bug in harker plots
Allow lists as inputs to most function key word arguments (e.g., P_bar etc.)
Added Fe redox to labeling
Add "Phase List" as a column header in outputs
Ensure consistent coloration when making multiple plots of phase abundance.

Version 0.3.2 (December 2025)
====================================
Updated installation functions for more streamlined installation process.

Version 0.3.0 (December 2025)
====================================
Update to the documentation of all functions. 
Units now included for all outputs and consistency in outputs between alphaMELTS and MAGEMin calculations.
Update MAGEMinCalc to v0.5.2 - use  `ptt.update_MAGEMinCalc()` function.

Version 0.2.55 (November 2025)
====================================
Outputs between MELTS and MAGEMin calculations are now consistent and indicated with correct units.
Update MAGEMinCalc to v0.5.1 - use `ptt.update_MAGEMinCalc()` function

Version 0.2.49 (October 2025)
====================================
In all notebooks ptt must now be imported as `import petthermotools`, `import PetThermoTools` is depreciated.
Major update to MAGEMinCalc - please use v0.5.0 and MAGEMin_C v2.0.6.
Phase Diagram code has better optimization and should have less calculation time
More plotting options added.

Version 0.2.40 (August 2025)
====================================
Adding melting calculations to MAGEMinCalc extension - works with MAGEMin_C up to v1.9.6.

Version 0.2.39 (July 2025)
====================================
Fix to minor bugs influencing crystallization calculations. Changes to the name of functions used for barometry calculations.
Please ensure you're updated to MAGEMinCalc v0.4.3 (see instructions) if running MAGEMin calculations.

Version 0.2.32 (March 28th 2025)
====================================
Changes to the underlying code and associated MAGEMinCalc library to improve the ease of integration between the packages. Adiabatic decompression melting in MAGEMinCalc is temporarily turn off.
MAGEMinCalc v0.4.1 https://github.com/gleesonm1/MAGEMinCalc/archive/refs/tags/v0.4.1.zip 

Version 0.2.30 (January 16th 2025)
====================================
Updates to MAGEMinCalc that allow greater integration with PetThermoTools and reduce the dependency on the pyMAGEMINcalc package.
Fixes to minor bugs included output of single equilibrium calculations.
MAGEMin_C v1.6.9 https://github.com/ComputationalThermodynamics/MAGEMin_C.jl/archive/refs/tags/v1.6.9.zip
MAGEMinCalc v0.3.8 https://github.com/gleesonm1/MAGEMinCalc/archive/refs/tags/v0.3.8.zip 

Version 0.2.25 (November 28th 2024)
=====================================
New 'supCalc' function added to work as a supplemental calculator
Automatically adds alphaMELTS for Python files to path if working in VICTOR.
Fixed other minor bugs.
pyMAGEMINcalc v0.0.8 https://github.com/gleesonm1/pyMAGEMINcalc/archive/refs/tags/v0.0.8.zip
MAGEMinCalc v0.3.4 https://github.com/gleesonm1/MAGEMinCalc/archive/refs/tags/v0.3.4.zip 
MAGEMin_C v1.3.5 (tested up to v1.5.9) https://github.com/ComputationalThermodynamics/MAGEMin_C.jl/archive/refs/tags/v1.3.5.zip 

Version 0.2.20 (October 20th 2024)
=====================================
Update to equilibrate_multi function to allow calculations with Holland et al. (2018) thermodynamic model.
Input code on ptt import to determine whether alphaMELTS files are in the Python path.
pyMAGEMINcalc v0.0.8 https://github.com/gleesonm1/pyMAGEMINcalc/archive/refs/tags/v0.0.8.zip
MAGEMinCalc v0.3.4 https://github.com/gleesonm1/MAGEMinCalc/archive/refs/tags/v0.3.4.zip 
MAGEMin_C v1.3.5 (tested up to v1.5.0) https://github.com/ComputationalThermodynamics/MAGEMin_C.jl/archive/refs/tags/v1.3.5.zip 

Version 0.2.18 (October 10th 2024)
=====================================
Updated installation guide for PetThermoTools with admin issues noted.
Update to equilibrate_multi to allow users to 'turn off' solid phases as requested.
Other minor bugs fixed.

Version 0.2.13 (September 17th 2024)
========================================
New installation guide for PetThermoTools and alphaMELTS for Python added.
equilibrate_multi now only returns one variable - a DataFrame containing all information about the calculation including the affinity of undersaturated phases.
A new approach to solving volatile saturation pressures has been added but has not been thoroughly tested.
Release for alphaMELTS workshop at GSA connects 2024.
pyMAGEMINcalc v0.0.7 https://github.com/gleesonm1/pyMAGEMINcalc/archive/refs/tags/v0.0.7.zip
MAGEMinCalc v0.3.2 https://github.com/gleesonm1/MAGEMinCalc/archive/refs/tags/v0.3.2.zip 
MAGEMin_C v1.3.5 (tested up to v1.5.0) https://github.com/ComputationalThermodynamics/MAGEMin_C.jl/archive/refs/tags/v1.3.5.zip 

Version 0.2.5 (January 11th, 2024)
========================================
Liquidus code updated to accept single (integer or float) inputs for initial temperature guess.
Path/crysatllization and degassing codes now return cumulate sum of fractionated phases - e.g., Results['Mass']['olivine1_cumsum']

Version 0.2.4 (January 9th, 2024)
======================================
Update to crystallization and decompression codes so that the output DataFrames are always the same length.

Version 0.2.1 (November 22nd, 2023)
=====================================
Repository name updated to PetThermoTools (Petrology Thermodynamics Tools).
Change-log and available functions linked to ReadTheDocs.
pyMAGEMINcalc, MAGEMinCalc, and MAGEMin_C versions are as listed above.

Version 0.1.29 (November 10th, 2023)
=======================================
Updated function name for mineral barometry code ('find_mineral_cosaturation').
Added schematic showing current available functions to ReadTheDocs.

For full functionality when using MAGEMin calculations, make sure PetThermoTools v0.1.29 is installed:
pyMAGEMINcalc v0.0.7 https://github.com/gleesonm1/pyMAGEMINcalc/archive/refs/tags/v0.0.7.zip
MAGEMinCalc v0.2.22 https://github.com/gleesonm1/MAGEMinCalc/archive/refs/tags/v0.2.22.zip 
MAGEMin_C v1.3.5 https://github.com/ComputationalThermodynamics/MAGEMin_C.jl/archive/refs/tags/v1.3.5.zip 

Version 0.1.28 (October 11th, 2023)
=======================================
Added adiabatic decompression melting function

Version 0.1.25 (August 8th, 2023)
======================================
Added both "comp" and "bulk" keyword arguments to multi_path for consistency with prior versions.