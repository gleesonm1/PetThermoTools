Naming conventions
==================

This page describes the naming terminology used in the output dataframes.

Thermodynamic properties
------------------------

The keys ``_prop`` (like ``liquid1_prop``) return the thermodynamic properties of each phase.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Name
     - Description
   * - ``mass_g_{phase}``
     - Mass in g of each phase (see equation below)
   * - ``rho_kg/m3_{phase}``
     - Density of phase in kg/m3
   * - ``v_cm3_{phase}``
     - Volume of phase in cm3
   * - ``G_J_{phase}``
     - Gibbs free energy in J of the phase
   * - ``H_J_{phase}``
     - Enthalpy in J of the phase
   * - ``S_J/K_{phase}``
     - Entropy in J/K of the phase
   * - ``Cp_J/kg/K_{phase}``
     - Specific heat capacity of the phase, J/(kg.K)
   * - ``dCpdT_J/(kg.K^2)_{phase}``
     - Partial derivative of specific heat capacity with respect to temperature, J/(kg.K^2)
   * - ``dVdT_cm^3/K_{phase}``
     - Partial derivative of phase volume with respect to temperature at constant pressure (cm^3/K)
   * - ``dPdT_bar/K_{phase}``
     - Partial derivative of pressure with respect to temperature at constant volume (bar K⁻¹)
   * - ``d2VdT2_cm^3/K^2_{phase}``
     - Second partial derivative of phase volume with respect to temperature at constant pressure (cm^3 K^-2)
   * - ``d2VdTdP_cm^3/(bar.K)_{phase}``
     - Mixed second partial derivative of phase volume with respect to temperature and pressure (cm^3 bar^-1 K^-1)
   * - ``d2VdP2_cm^3/bar^2_{phase}``
     - Second partial derivative of phase volume with respect to pressure at constant temperature (cm^3 bar^-2)

Then per oxide and per phase, the following are given 


.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Name
     - Description
   * - ``activity_{oxide)_{phase}``
     - Activity of each oxide in the phase (e.g. activity_Na2O_Liq is activity of Na2O in liquid)
  * - ``activity0_{oxide)_{phase}``
     - Standard State activity of each oxide in the phase (e.g. activity0_Na2O_Liq is standard state activity of Na2O in liquid)
  * - ``mu_{oxide)_{phase}``
     - Chemical potential of the oxide in the phase, Na2O in liquid, in J/mol)
  * - ``mu0_{oxide)_{phase}``
     - Standard State Chemical potential of the oxide in the phase, Na2O in liquid, in J/mol)
  * - ``X_{oxide)_{phase}``
     - Mole fraction of the oxide in the phase. 
