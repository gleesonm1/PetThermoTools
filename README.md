## pyMELTScalc

Petrological calculations using the MELTS or Holland (MAGEMin) thermodynamic databases.

Before this package can be used, the Julia code used to compile and run the MAGEMin calculations must be installed. To do this:
1. Install Julia.
2. open Julia and type
a. import Pkg
b. using Pkg
c. Pkg.add(url = "https://github.com/gleesonm1/MAGEMinCalc")
d. ]
e. precompile

