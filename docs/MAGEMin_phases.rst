#################################
Mineral phases available in MAGEMin calculations
#################################

A couple of common questions we receive are "Can MELTS/MAGEMin predict X phase?" and "Can MELTS/MAGEMin model Y component in X phase?"

To try and help all users understand the capabilities and limitations of each thermodynamic, the available phases of each model used in MAGEMin (and available through PetThermoTools) are listed below. 

.. note:: 
   This list is adapted from the the MAGEMin documentation https://computationalthermodynamics.github.io/MAGEMin_C.jl/dev/



Green et al. (2025) - named `ig` in MAGEMin
=====================

:math:`K_2O - Na_2O - CaO - FeO - MgO - Al_2O_3 - SiO_2 - H_2O - TiO_2 - O - Cr_2O_3` chemical system

.. note::
    In PetThermoTools we specify an fO2 buffer and/or an Fe redox state so the `O` content of the system is defined indirectly by one of these two parameters.

Pure phases
---------------------
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
----------------------

.. note::
    All Fe is 2+ valence unless specified. End member information taken from https://hpxeosandthermocalc.org/the-hpx-eos/the-hpx-eos-families/hpx-eos-igneous-sets

* **olivine**
    * mont :math:`CaMgSiO_4`
    * fa :math:`Fe_2SiO_4`
    * fo :math:`Mg_2SiO_4`
    * cfm :math:`MgFeSiO_4` (ordered - Mg in M1, Fe in M2)
* **feldspar**
    * ab :math:`NaAlSi_3O_8`
    * san :math:`KAlSi_3O_8`
    * an :math:`CaAl_2Si_2O_8`
* **clinopyroxene**
    * di :math:`CaMgSi_2O_6`
    * cfs :math:`Fe_2Si_2O_6`
    * cats :math:`CaAl_2SiO_6`
    * crdi :math:`CaCrAlSiO_6`
    * cess :math:`CaFe^{3+}AlSiO_6`
    * cbuf :math:`Ca_{1.5}AlTi_{0.5}SiO_6`
    * jd :math:`NaAlSi_2O_6`
    * cen :math:`Mg_2Si_2O_6`
    * cfm :math:`MgFeSi_2O_6` (ordered - Mg in M1, Fe in M2)
    * kjd :math:`KAlSi_2O_6`
* **orthopyroxene**
    * en :math:`Mg_2Si_2O_6`
    * fs :math:`Fe_2Si_2O_6`
    * fm :math:`MgFeSi_2O_6` (ordered - Mg in M1, Fe in M2)
    * odi :math:`CaMgSi_2O_6`
    * mgts :math:`MgAl_2SiO_6`
    * cren :math:`MgCrAlSiO_6`
    * obuf :math:`Mg_{1.5}AlTi_{0.5}SiO_6`
    * mess :math:`MgFe^{3+}AlSiO_6`
    * ojd :math:`NaAlSi_2O_6`
* **garnet**
    * py :math:`Mg_3Al_2Si_3O_{12}`
    * alm :math:`Fe_3Al_2Si_3O_{12}`
    * gr :math:`Ca_3Al_2Si_3O_{12}`
    * andr :math:`Ca_3Fe^{3+}_2Si_3O_{12}`
    * knor :math:`Mg_3Cr_2Si_3O_{12}`
    * tig :math:`Mg_{3.5}AlTi_{0.5}Si_3O_{12}`
* **clino-amphibole**
    * tr :math:`Ca_2Mg_5Si_8O_{22}(OH)_2`
    * tsm :math:`Ca_2Mg_3Al_4Si_6O_{22}(OH)_2`
    * prgm :math:`NaCa_2Mg_4Al_3Si_6O_{22}(OH)_2`
    * cumm :math:`Mg_7Si_8O_{22}(OH)_2`
    * grnm :math:`Fe_7Si_8O_{22}(OH)_2`
    * a :math:`Mg_3Fe_4Si_8O_{22}(OH)_2` (ordered - Fe in M4 and M2, Mg in M1 and M3)
    * b :math:`Mg_2Fe_5Si_8O_{22}(OH)_2` (ordered - Fe in M1, M3, and M4, Mg in M2)
    * mrb :math:`Na_2Mg_3Fe^{3+}_2Si_8O_{22}(OH)_2`
    * kprg :math:`KCa_2Mg_4Al_3Si_6O_{22}(OH)_2`
    * tts :math:`Ca_2Mg_3Al_2Ti_2Si_6O_{24}`
* **spinel**
    * nsp :math:`MgAl_2O_4` (ordered)
    * isp :math:`MgAl_2O_4` (ordered, inverse)
    * nhc :math:`FeAl_2O_4` (ordered)
    * ihc :math:`FeAl_2O_4` (ordered, inverse)
    * nmt :math:`FeFe^{3+}_2O_4` (ordered)
    * imt :math:`FeFe^{3+}_2O_4` (ordered, inverse)
    * picr :math:`MgCr_2O_4`
    * usp :math:`Fe_2TiO_4`
* **biotite**
    * phl :math:`KMg_3AlSi_3O_{10}(OH)_2`
    * annm :math:`KFe_3AlSi_3O_{10}(OH)_2`
    * obi :math:`KMg_2FeAlSi_3O_{10}(OH)_2` (ordered - Fe in M3 site)
    * est :math:`KMg_2Al_3Si_2O_{10}(OH)_2`
    * tbi :math:`KMg_2AlSi_3TiO_{12}`
    * fbi :math:`KMg_2Al_2Fe^{3+}Si_2O_{10}(OH)_2`
* **cordierite**
    * crd :math:`Mg_2Al_4Si_5O_{18}`
    * fcrd :math:`Fe_2Al_4Si_5O_{18}`
    * hcrd :math:`Mg_2Al_4Si_5O_{17}(OH)_2`
* **epidote**
    * cz :math:`Ca_2Al_3Si_3O_{12}(OH)`
    * ep :math:`Ca_2FeAl_2Si_3O_{12}(OH)`
    * fep :math:`Ca_2Fe_2AlSi_3O_{12}(OH)`
* **ilmenite**
    * oilm :math:`FeTiO_3` (Fe in A site, Ti in B)
    * dilm :math:`FeTiO_3` (Fe and Ti split across A and B)
    * hem :math:`Fe^{3+}_2O_3`
    * ogk :math:`MgTiO_3` (Mg in A site, Ti in B)
    * dgk :math:`MgTiO_3` (Mg and Ti split across A and B)
* **muscovite**
    * mu :math:`KAl_3Si_3O_{12}(OH)_2`
    * cel :math:`KMgAlSi_4O_{10}(OH)_2`
    * fcel :math:`KFeAlSi_4O_{10}(OH)_2`
    * pa :math:`NaAl_3Si_3O_{10}(OH)_2`
    * mat :math:`CaAl_4Si_2O_{10}(OH)_2`
    * fmu :math:`KAl_2Fe^{3+}Si_3O_{12}(OH)2`
* **silicate melt**
    * q4L :math:`Si_4O_8`
    * slL :math:`Al_2SiO_5`
    * wo1L :math:`CaSiO3`
    * fo2L :math:`Mg_4Si_2O_8`
    * fa2L :math:`Fe_4Si_2O_8`
    * jdL :math:`NaAlSi_2O_6`
    * hml :math:`Fe^{3+}O_{1.5}`
    * ekl :math:`CrO_{1.5}`
    * tiL :math:`TiO_2`
    * kjL :math:`KAlSi_2O_6`
    * ctL :math:`CaAl_2SiO_6`
    * h2O1L :math:`H_2O`
* **fluid**


Weller et al. (2024) - named `igad` in MAGEMin
=====================

:math:`K_2O - Na_2O - CaO - FeO - MgO - Al_2O_3 - SiO_2 - TiO_2 - O - Cr_2O_3` chemical system

.. note::
    In PetThermoTools we specify an fO2 buffer and/or an Fe redox state so the `O` content of the system is defined indirectly by one of these two parameters. Additionally, the Weller et al. (2024) model is currently anhydrous.

Pure phases
--------------------
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
--------------------

.. note::
    All Fe is 2+ valence unless specified. End member information taken from https://hpxeosandthermocalc.org/the-hpx-eos/the-hpx-eos-families/hpx-eos-igneous-sets

* **olivine**
    * mont :math:`CaMgSiO_4`
    * fa :math:`Fe_2SiO_4`
    * fo :math:`Mg_2SiO_4`
    * cfm :math:`MgFeSiO_4` (ordered - Mg in M1, Fe in M2)
* **feldspar**
    * ab :math:`NaAlSi_3O_8`
    * san :math:`KAlSi_3O_8`
    * an :math:`CaAl_2Si_2O_8`
* **clinopyroxene**
    * di :math:`CaMgSi_2O_6`
    * cfs :math:`Fe_2Si_2O_6`
    * cats :math:`CaAl_2SiO_6`
    * crdi :math:`CaCrAlSiO_6`
    * cess :math:`CaFe^{3+}AlSiO_6`
    * cbuf :math:`Ca_{1.5}AlTi_{0.5}SiO_6`
    * jd :math:`NaAlSi_2O_6`
    * cen :math:`Mg_2Si_2O_6`
    * cfm :math:`MgFeSi_2O_6` (ordered - Mg in M1, Fe in M2)
    * kjd :math:`KAlSi_2O_6`
* **orthopyroxene**
    * en :math:`Mg_2Si_2O_6`
    * fs :math:`Fe_2Si_2O_6`
    * fm :math:`MgFeSi_2O_6` (ordered - Mg in M1, Fe in M2)
    * odi :math:`CaMgSi_2O_6`
    * mgts :math:`MgAl_2SiO_6`
    * cren :math:`MgCrAlSiO_6`
    * obuf :math:`Mg_{1.5}AlTi_{0.5}SiO_6`
    * mess :math:`MgFe^{3+}AlSiO_6`
    * ojd :math:`NaAlSi_2O_6`
* **garnet**
    * py :math:`Mg_3Al_2Si_3O_{12}`
    * alm :math:`Fe_3Al_2Si_3O_{12}`
    * gr :math:`Ca_3Al_2Si_3O_{12}`
    * andr :math:`Ca_3Fe^{3+}_2Si_3O_{12}`
    * knor :math:`Mg_3Cr_2Si_3O_{12}`
    * tig :math:`Mg_{3.5}AlTi_{0.5}Si_3O_{12}`
* **spinel**
    * nsp :math:`MgAl_2O_4` (ordered)
    * isp :math:`MgAl_2O_4` (ordered, inverse)
    * nhc :math:`FeAl_2O_4` (ordered)
    * ihc :math:`FeAl_2O_4` (ordered, inverse)
    * nmt :math:`FeFe^{3+}_2O_4` (ordered)
    * imt :math:`FeFe^{3+}_2O_4` (ordered, inverse)
    * picr :math:`MgCr_2O_4`
    * usp :math:`Fe_2TiO_4`
* **ilmenite**
    * oilm :math:`FeTiO_3` (Fe in A site, Ti in B)
    * dilm :math:`FeTiO_3` (Fe and Ti split across A and B)
    * hem :math:`Fe^{3+}_2O_3`
    * ogk :math:`MgTiO_3` (Mg in A site, Ti in B)
    * dgk :math:`MgTiO_3` (Mg and Ti split across A and B)
* **nepheline**
    * neN :math:`Na_4Al_4Si_4O_{16}`
    * neS :math:`Na_3Al_3Si_5O_{16}` (vacancy in A2 site)
    * neK :math:`K_4Al_4Si_4O_{16}`
    * neO :math:`Na_3KAl_4Si_4O_{16}` (ordered - Na in A1, K in A2)
    * neC :math:`Na_2CaAl_4Si_4O_{16}` (vacancy in A2 site)
    * neF :math:`Na_4Fe^{3+}_4Si_4O_{16}`
* **kalsilite**
    * nks :math:`NaAlSiO_4`
    * kls :math:`KAlSiO_4`
* **leucite**
    * nlc :math:`NaAlSi_2O_6`
    * lc :math:`KAlSi_2O_6`
* **melilite**
    * geh :math:`Ca_2Al_2SiO_7`
    * ak :math:`Ca_2MgSi_2O_7`
    * fak :math:`Ca_2FeSi_2O_7`
    * nml :math:`NaCaAl_2SiO_7`
    * fge :math:`Ca_2Fe^{3+}AlSiO_7`
* **silicate melt**
    * q3L :math:`Si_3O_6`
    * sl1L :math:`Al_2SiO_5`
    * wo1L :math:`CaSiO_3`
    * fo2L :math:`Mg_4Si_2O_8`
    * fa2L :math:`Fe_4Si_2O_8`
    * nml :math:`NaSi_{0.5}O_{1.5}`
    * hml :math:`Fe^{3+}O_{1.5}`
    * ekl :math:`CrO_{1.5}`
    * tiL :math:`TiO_2`
    * kmL :math:`KSi_{0.5}O_{1.5}`
    * anL :math:`CaAl_2Si_2O_8` (ordered)
    * ab1L :math:`NaAlSi_3O_8` (ordered)
    * kfL :math:`KAlSi_3O_8` (ordered)
    * enL :math:`Mg_2Si_2O_6` (ordered)
