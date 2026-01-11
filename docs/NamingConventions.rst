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
   * - ``g_J_{phase}``
     - Gibbs free energy in J of the phase
   * - ``h_J_{phase}``
     - Enthalpy in J of the phase
   * - ``s_J/K_{phase}``
     - Entropy in J/K of the phase
   * - ``cp_J/kg/K_{phase}``
     - Specific heat capacity of the phase (J kg-1 K-1)
   * - ``dcpdt_J/kg/K_{phase}``
     - Partial derivative of specific heat capacity with respect to temperature (J kg⁻¹ K⁻¹)
   * - ``dvdt_cm3/K_{phase}``
     - Partial derivative of phase volume with respect to temperature at constant pressure (cm³ K⁻¹)
   * - ``dpdt_bar/K_{phase}``
     - Partial derivative of pressure with respect to temperature at constant volume (bar K⁻¹)
