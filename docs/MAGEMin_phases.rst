#################################
Mineral phases available in MAGEMin calculations
#################################

A couple of common questions we receive are "Can MELTS/MAGEMin predict X phase?" and "Can MELTS/MAGEMin model Y component in X phase?"

To try and help all users understand the capabilities and limitations of each thermodynamic, the available phases of each model used in MAGEMin (and available through PetThermoTools) are listed below. 

.. note:: 
   This list is adapted from the the MAGEMin documentation https://computationalthermodynamics.github.io/MAGEMin_C.jl/dev/


====================
Green et al. (2025) - named `ig` in MAGEMin
=====================

:math:`K_2O - Na_2O - CaO - FeO - MgO - Al_2O_3 - SiO_2 - H_2O - TiO_2 - O - Cr_2O_3` chemical system

.. note::
    In PetThermoTools we specify an fO2 buffer and/or an Fe redox state so the `O` content of the system is defined indirectly by one of these two parameters.

Pure phases
=====================
* **quartz** :math:`SiO_2`
* **cristobalite** :math:`SiO_2`
* **tridymite** :math:`SiO_2`
* **coesite** :math:`SiO_2`
* **stishovite** :math:`SiO_2`
* **kyanite** :math:`Al_2O_5`
* **sillimanite** :math:`Al_2O_5`
* **andalusite** :math:`Al_2O_5`
* **rutile** :math:`TiO_2`
* **sphene** :math:`CaTiSiO_5`

Solution phases
======================
* **olivine**
* **feldspar**
* **clinopyroxene**
* **orthopyroxene**
* **garnet**
* **clino-amphibole**
* **spinel**
* **biotite**
* **cordierite**
* **epidote**
* **ilmenite**
* **muscovite**
* **fluid**
* **silicate melt**

====================
Weller et al. (2024) - named `igad` in MAGEMin
=====================

:math:`K_2O - Na_2O - CaO - FeO - MgO - Al_2O_3 - SiO_2 - TiO_2 - O - Cr_2O_3` chemical system

.. note::
    In PetThermoTools we specify an fO2 buffer and/or an Fe redox state so the `O` content of the system is defined indirectly by one of these two parameters. Additionally, the Weller et al. (2024) model is currently anhydrous.

Pure phases
=====================
* **quartz** :math:`SiO_2`
* **cristobalite** :math:`SiO_2`
* **tridymite** :math:`SiO_2`
* **coesite** :math:`SiO_2`
* **stishovite** :math:`SiO_2`
* **kyanite** :math:`Al_2O_5`
* **sillimanite** :math:`Al_2O_5`
* **andalusite** :math:`Al_2O_5`
* **rutile** :math:`TiO_2`
* **sphene** :math:`CaTiSiO_5`

Solution phases
======================
* **olivine**
* **feldspar**
* **clinopyroxene**
* **orthopyroxene**
* **garnet**
* **spinel**
* **ilmenite**
* **nepheline**
* **kalsilite**
* **leucite**
* **melilite**
* **silicate melt**
