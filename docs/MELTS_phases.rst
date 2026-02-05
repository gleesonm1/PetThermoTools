#################################
Mineral phases available in MELTS
#################################

A couple of common questions we receive are "Can MELTS predict X phase?" and "Can MELTS model Y component in X phase?"

To try and help all users understand the capabilities and limitations of MELTS, the available phases (and their end-members for solid solution phases) are listed below. 

.. note:: 
   This list is adapted from the `original MagmaSource forum <https://magmasource.caltech.edu/forum/index.php?topic=85.0>`_.

Solid solution phases
=====================

* **olivine**
    * tephroite: :math:`Mn_2SiO_4`
    * fayalite: :math:`Fe_2SiO_4`
    * Co-olivine: :math:`Co_2SiO_4`
    * Ni-olivine: :math:`Ni_2SiO_4`
    * monticellite: :math:`CaMgSiO_4`
    * forsterite: :math:`Mg_2SiO_4`

* **garnet**
    * almandine :math:`Fe^{2+}_3Al_2Si_3O_{12}`
    * grossular :math:`Ca_3Al_2Si_3O_{12}`
    * pyrope :math:`Mg_3Al_2Si_3O_{12}`

* **melilite**
    * akermanite :math:`Ca_2MgSi_2O_7`
    * gehlenite :math:`Ca_2Al_2SiO_7`
    * Fe-akermanite :math:`Ca_2Fe^{2+}Si_2O_7`
    * soda-melilite :math:`Na_2Si_3O_7`

* **orthopyroxene**
    * diopside :math:`CaMgSi_2O_6`
    * clinoenstatite :math:`Mg_2Si_2O_6`
    * hedenbergite :math:`CaFe^{2+}Si_2O_6`
    * alumino-buffonite :math:`CaTi_{0.5}Mg_{0.5}AlSiO_6`
    * buffonite :math:`CaTi_{0.5}Mg_{0.5}Fe^{3+}SiO_6`
    * essenite :math:`CaFe^{3+}AlSiO_6`
    * jadeite :math:`NaAlSi_2O_6`

* **clinopyroxene**
    * diopside :math:`CaMgSi_2O_6`
    * clinoenstatite :math:`Mg_2Si_2O_6`
    * hedenbergite :math:`CaFe^{2+}Si_2O_6`
    * alumino-buffonite :math:`CaTi_{0.5}Mg_{0.5}AlSiO_6`
    * buffonite :math:`CaTi_{0.5}Mg_{0.5}Fe^{3+}SiO_6`
    * essenite :math:`CaFe^{3+}AlSiO_6`
    * jadeite :math:`NaAlSi_2O_6`

.. note::
    The compositions of the essenite, buffonite and alumino-buffonite end members mean that calculations for bulk compositions that contain just two of Fe2O3, TiO2 and Al2O3 are likely to fail (e.g. no TiO2 in the bulk requires that buffonite and alumino-buffonite cancel out exactly). It is safer to add trace amounts of these components to your calculations.

* **orthoamphibole**
    * cummingtonite :math:`Mg_7Si_8O_{22}(OH)_2`
    * grunerite :math:`Fe^{2+}_7Si_8O_{22}(OH)_2`
    * tremolite :math:`Ca_2Mg_5Si_8O_{22}(OH)_2`

* **clinoamphibole**
    * cummingtonite :math:`Mg_7Si_8O_{22}(OH)_2`
    * grunerite :math:`Fe^{2+}_7Si_8O_{22}(OH)_2`
    * tremolite :math:`Ca_2Mg_5Si_8O_{22}(OH)_2`

* **hornblende**
    * pargasite :math:`NaCa_2Mg_4AlAl_2Si_6O_{22}(OH)_2`
    * ferropargasite :math:`NaCa_2Fe^{2+}_4AlAl_2Si_6O_{22}(OH)_2`
    * magnesiohastingsite :math:`NaCa_2Mg_4Fe^{3+}Al_2Si_6O_{22}(OH)_2`

.. note::
    In some MELTS software there is just one amphibole solid solution rather than separate clino- and ortyho- phases.

* **biotite**
    * annite :math:`KFe^{2+}_3Si_3AlO_{10}(OH)_2`
    * phlogopite :math:`KMg_3Si_3AlO_{10}(OH)_2`

* **feldspar**
    * albite :math:`NaAlSi_3O_8`
    * anorthite :math:`CaAl_2Si_2O_8`
    * sanidine :math:`KAlSi_3O_8`

* **nepheline**
    * Na-nepheline :math:`Na_4Al_4Si_4O_{16}`
    * K-nepheline :math:`K_4Al_4Si_4O_{16}`
    * vc-nepheline :math:`Na_3Al_3Si_5O_{16}`
    * Ca-nepheline :math:`CaNa_2Al_4Si_4O_{16}`

* **leucite**
    * leucite :math:`KAlSi_2O_6`
    * analcime :math:`NaAlSi_2O_5(OH)_2`
    * na-leucite :math:`NaAlSi_2_O_6`

* **spinel**
    * chromite :math:`Fe^{2+}Cr_2O_4`
    * hercynite :math:`Fe^{2+}Al_2O_4`
    * magnetite :math:`Fe^{2+}Fe^{3+}_2O_4`
    * spinel :math:`MgAl_2O_4`
    * ulvospinel :math:`Fe^{2+}_2TiO_4`

* **rhm-oxide**
    * giekelite :math:`MgTiO_3`
    * hematite :math:`Fe^{3+}O_3`
    * ilmenite :math:`Fe^{2+}TiO_3`
    * pyrophanite :math:`MnTiO_3`

* **ortho-oxide**
    * psuedobrookite :math:`Fe^{3+}_2TiO_5`
    * ferropseudobrookite :math:`Fe^{2+}Ti_2O_5`
    * karrooite :math:`MgTi_2O_5`

* **alloy-solid**
    * Fe-metal :math:`Fe`
    * Ni-metal :math:`Ni`

* **alloy-liquid**
    * Fe-liquid :math:`Fe`
    * Ni-liquid :math:`Ni`


Pure phases included in all versions of MELTS
=====================
* **sphene** :math:`CaTiSiO_5`
* **aenigmatite** :math:`Na_2Fe_5TiSi_6_O_{20}`
* **muscovite** :math:`KAl_2Si_3AlO_{10}(OH)_2`
* **quartz** :math:`SiO_2`
* **tridymite** :math:`SiO_2`
* **cristobalite** :math:`SiO_2`
* **corundum** :math:`Al_2O_3`
* **rutile** :math:`TiO_2`
* **perovskite** :math:`CaTiO_3`
* **whitlockite** :math:`Ca_3(PO_4)_2`
* **apatite** :math:`Ca_4(PO_4)_3(OH)`


Additional pure phases included in alphaMELTS
=====================
* **andalusite** :math:`Al_2SiO_5`
* **sillimanite** :math:`Al_2SiO_5`
* **kyanite** :math:`Al_2SiO_5`
* **coesite** :math:`SiO_2`
* **talc** :math:`Mg_3Si_4O_{10}(OH)_2`
* **anthophyllite** :math:`Mg_7Si_8_O_{22}(OH)_2`
* **cordierite** :math:`Mg_2Al_4Si_5O_{18}`
* **lawsonite** :math:`CaAl_2Si_2_O_6(OH)_4`
* **antigorite** :math:`Mg_{48}Si_{34}O_{85}(OH)_{62}`
* **chrysotile** :math:`Mg_{3}Si_{2}O_{5}(OH)_{4}`
* **brucite** :math:`Mg(OH)_2`