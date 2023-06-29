============
Installation
============

There are three steps to installing pyMELTScalc locally on your machine.

Step 1 - Install Python
======================

This may seem obvious, but the first step is to obtain Python3 (tested on V3.7). If you haven't used python before, we recomend installing it through anaconda.
 `anaconda3 <https://www.anaconda.com/products/individual>`_.

Step 2 - Download the alphaMELTS for MATLAB/Python source files
================================================================

pyMELTScalc itself doesn't run any thermodynamic calculations. Instead, pyMELTScalc simply organises your inputs to so that the calculations will run in alphaMELTS for Python without you having to write out a new scipt every time you want to do something new! It also allows you to automate a huge number of calculations (e.g., enabling Monte Carlo simulations).
Therefore, before any calculations are performed in pyMELTScalc users need to download the alphaMELTS for MATLAB/Python files from here: https://magmasource.caltech.edu/gitlist/MELTS_Matlab.git/ unzip the package and save it somewhere you'll remember!

Once this has been completed before any Python script using pyMELTScalc can be run, the the location of the alphaMELTS for MATLAB/Python files (found under the package subfolder in the above link) must be added to the Python path. It is possible to add a folder permanently to the Python path, but we recommend simply using the following code at the start of each Python script using pyMELTScalc:

.. code-block:: python

   import sys
   sys.path.append(r'LOCATION OF MELTS/PACKAGE files')


Step 3 - Install pyMELTScalc
============================

pyMELTScalc can be installed using pip in one line. If you are using a terminal, enter:

.. code-block:: python

   pip install pyMELTScalc

If you are using Jupyter Notebooks or Jupyter Lab, you can also install it by entering the following code into a notebook cell (note the !):

.. code-block:: python

   !pip install pyMELTScalc

You then need to import pyMELTScalc into the script you are running code in. In all the examples, we import pyMELTScalc as M:

.. code-block:: python

   import pyMELTScalc as M

This means any time you want to call a function from pyMELTScalc, you do M.function_name.



Updating
========

To upgrade to the most recent version of pyMELTScalc, type the following into terminal:

.. code-block:: python

   pip install pyMELTScalc --upgrade

Or in your Jupyter environment:

.. code-block:: python

   !pip install pyMELTScalc --upgrade


For maximum reproducability, you should state which version of pyMELTScalc you are using. If you have imported pyMELTScalc as M, you can find this using:

.. code-block:: python

    M.__version__