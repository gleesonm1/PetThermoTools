## pyMELTScalc

Petrological calculations using the MELTS or thermodynamic databases. Through integration with the pyMAGEMINcalc repository, calculations can also be performed using the thermodynamic database of Holland et al. (2018).

# Using MELTS
As a default, this package uses the MELTS algorithms and the alphaMELTS for Python package, developed by Dr Paula Antoshechkina, to perform the thermodynamic calculations. Therefore, it is neccessary for any user to first download the alphaMELTS for Python package. Luckily, version 1.2.0beta of MELTS for Python is included in this repository and we have extensively checked that this version of MELTS for Python works with out code. Future versions of MELTS for Python might be available via:
https://magmasource.caltech.edu/gitlist/MELTS_Matlab.git/
although we caution that newer versions of MELTS for Python may not be compatible with pyMELTScalc.

Prior to any calculations, the user must add the location of the meltsdynamic.py file to the path using import sys; sys.path.append(dir).

# Using MAGEMin
This package also provides the user to perform calculations using the thermodynamic database of Holland et al. (2018). This is made possible through the recent development and release of MAGEMin (https://github.com/ComputationalThermodynamics; https://doi.org/10.1029/2022GC010427), with Julia functions called from Python used to run the calculations.

Before this package can be used, the Julia code used to compile and run the MAGEMin calculations, which is hosted in a separate repository, must be imported and added to Julia.
First, you will need to open Python and install the PyJulia packages by running the following lines:
import julia
julia.install()

Then, the MAGEMinCalc functions can be installed via:

1. Install Julia (https://julialang.org/downloads/).

2. Open Python and install Julia via:
a."pip install Julia"
b. import julia
c. julia.install()

3. Run Julia and add the MAGEMinCalc package via the following commands:
a. import Pkg
b. using Pkg
c. Pkg.add(url = "https://github.com/gleesonm1/MAGEMinCalc")

Finally, the pyMAGEMINcalc package can be installed in a Jupyter terminal via:
!pip install "https://github.com/gleesonm1/pyMAGEMINcalc/archive/refs/tags/v0.0.1.zip"

# Current calculations
At the moment, only a small selection of calculations are possible in pyMELTScalc, but this is expected to expand rapidly in the future.

Currently, users can:
1. Search for melt liquidus temperatures.
2. Run isobaric, polybaric, and isochoric crystallisation models.
3. Perform barometry (and hygrometry) calculations by searching for the P, H2O and fO2 values required for multi-phase saturation to occur.
