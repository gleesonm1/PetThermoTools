============
Using MAGEMin
=============
Before any calculations can be performed using the MAGEMin software you must have Python installed on your computed and PetThermoTools installed via pip.

- First ensure that they have a local Julia installation (https://julialang.org/downloads/).
- Next, users should run Julia and enter the following commands:

- import Pkg 
- using Pkg 
- Pkg.add(url = "https://github.com/gleesonm1/MAGEMinCalc")

- Once this has been done, users need to open Python and install Julia within the Python environment via pip. This can be done in the terminal:

.. code-block:: python

   pip install Julia

or within a Jupyter environment:

.. code-block:: python

   !pip install Julia

- Following this, users should run the following code in Python:

.. code-block:: python

   import julia 
   julia.install()

- At this point the installation is nearly complete, users simply need to install the pyMAGEMINcalc package in Python via the terminal:

.. code-block:: python

   pip install "https://github.com/gleesonm1/pyMAGEMINcalc/archive/refs/tags/v0.0.7.zip"

or a Jupyter environment:

.. code-block:: python

   !pip install "https://github.com/gleesonm1/pyMAGEMINcalc/archive/refs/tags/v0.0.7.zip"

