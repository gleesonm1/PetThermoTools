========================
Installation & Updating
========================

Installation
============

First, obtain Python3 (tested on V3.7). If you haven't used python before, we recomend installing it through anaconda.
 `anaconda3 <https://www.anaconda.com/products/individual>`_.

All of the You Tube examples shown here use Jupyter Lab, which is an easy-friendly editor that will be installed along with Anaconda.

pyMELTScalc can be installed using pip in one line. If you are using a terminal, enter:

.. code-block:: python

   pip install pyMELTScalc

If you are using Jupyter Notebooks or Jupyter Lab, you can also install it by entering the following code into a notebook cell (note the !):

.. code-block:: python

   !pip install pyMELTScalc

You then need to import pyMELTScalc into the script you are running code in. In all the examples, we import pyMELTScalc as ss:

.. code-block:: python

   import pyMELTScalc as ss

This means any time you want to call a function from pyMELTScalc, you do ss.function_name.



Updating
========

To upgrade to the most recent version of pyMELTScalc, type the following into terminal:

.. code-block:: python

   pip install pyMELTScalc --upgrade

Or in your Jupyter environment:

.. code-block:: python

   !pip install pyMELTScalc --upgrade


For maximum reproducability, you should state which version of pyMELTScalc you are using. If you have imported pyMELTScalc as pt, you can find this using:

.. code-block:: python

    ss.__version__