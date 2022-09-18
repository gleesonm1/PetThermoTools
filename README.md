## pyMELTScalc

Petrological calculations using the MELTS or Holland (MAGEMin) thermodynamic databases.

# Using MELTS
As a default, this package uses the MELTS algorithms and the alphaMELTS for Python package, developed by Dr Paula Antoshechkina, to perform the thermodynamic calculations. Therefore, it is neccessary for any user to first download the alphaMELTS for Python package from:
https://magmasource.caltech.edu/gitlist/MELTS_Matlab.git/

Prior to any calculations, the user must then add the location of the meltsdynamic.py file using import sys; sys.path.append(dir).

# Using MAGEMin
This package also provides the user to perform calculations using the thermodynamic database of Holland et al. (2018). This is made possible through the recent development and release of MAGEMin (https://github.com/ComputationalThermodynamics; https://doi.org/10.1029/2022GC010427), with Julia functions called from Python used to run the calculations.

Before this package can be used, the Julia code used to compile and run the MAGEMin calculations, which is hosted in a separate repository, must be imported and added to Julia. There are currently some minor problems that we are working through, but it should be possible to add the MAGEMinCalc package should by:
1. Install Julia (https://julialang.org/downloads/).
2. Run Julia and add the MAGEMinCalc package via the following commands:
a. import Pkg

b. using Pkg

c. Pkg.add(url = "https://github.com/gleesonm1/MAGEMinCalc")

You will then need to open Python and install the PyJulia packages by running the following lines:
import julia
julia.install()

# Current calculations
At the moment, only a small selection of calculations are possible in pyMELTScalc, but this is expected to expand rapidly in the future.

Currently, users can:
1. Search for melt liquidus temperatures.
2. Run isobaric, polybaric, and isochoric crystallisation models.
3. Perform barometry (and hygrometry) calculations by searching for the P, H2O and fO2 values required for multi-phase saturation to occur.
