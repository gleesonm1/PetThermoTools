#################
PetThermoTools in the Classroom!
#################

We have developed a number of materials for our classes to use PetThermoTools in intro, intermediate and advanced level classes. We outline these excersizes below and provide the required notebooks, student handouts, and answer keys. 

Investigating viscosity, liquidus temperatures, and fractional crystallization trends
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This practical was developed for EPS114/214 at UC Berkeley - a mixed undergraduate/graduate class in Petrology. In the three lectures prior, we discuss the liquidus and solidus, fractional crystallization, and the physical properties of magma like viscosity. In the practical before, we do the Skittles/M&M magma chamber so they have some intuition that removing a phase poor in element X drives it up in the remaining melt, removing a phase rich in element X drives it down in the remaining melt. This practical takes over from this. The learning objectives are:

-	Part 1: Understand what controls viscosity, and the relative effect of water, temperature, and composition 

-	Part 2: Understand the concept of a liquidus, and play around with parameters that control the liquidus phase. 

-	Part 3: Understand the influence of oxygen fugacity and water content on the development of tholeiitic vs. calc alkaline magmas. 

Students do not have to have any python experience, or even a local python installation. In the week prior to this class slot, I make them all sign up for a VICTOR account, and then they just have to download the notebook and upload it to VICTOR (https://victor.ldeo.columbia.edu/). I also give the students extra credit if they get PetThermoTools running locally. The first few questoins just require them to run cells and change inputs. We take about 40 minutes in class on this, and they complete the rest for homework (mostly just copying and pasting cells and editing code slightly). 

The Jupyter Notebook for students to work from is here:
https://github.com/gleesonm1/PetThermoTools/blob/master/docs/teaching_materials/FractionalCrystallization_Viscosity_CalcAlkalineTrends/PropertiesofMagmas_MELTSCrystallizationExamples_Final.ipynb

The homework exersize for students to fill in is here:
https://github.com/gleesonm1/PetThermoTools/blob/master/docs/teaching_materials/FractionalCrystallization_Viscosity_CalcAlkalineTrends/HW2-%20Fractional%20Crystallization_StudentVersion.docx

And the answer key here:
https://github.com/gleesonm1/PetThermoTools/blob/master/docs/teaching_materials/FractionalCrystallization_Viscosity_CalcAlkalineTrends/HW2-%20Fractional%20Crystallization_InstructorVersion.docx


Investigating Volatile Solubility and degassing behavoir
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This practical was developed for EPS114/214 at UC Berkeley - a mixed undergraduate/graduate class in Petrology. In the lectures before, we discuss how volatiles dissolve in magmas, discuss different solubility models, and the importance of volatiles as a driving force of eruptions. The aim of this practical was first use VESical (Iacovino et al. 2021) to get students to investigate how solubility models differ, and understand what non-ideality looks like (we discussed it in lectures). We use IaconoMarziano and Dixon as neither of these requires students to use a local thermoengine installation. 
Then students use PetThermoTools to investigate what happens to volatiles from the mantle to the surface. They start by ascending a magma to see when volatile saturation initiates. Then, they let the same magma undergo fractional crystallization, noting the differential behavoir of H2O and CO2. Finally, they ascend it towards the surface to see the large change in volatile volume (i.e. the driving force of volcanic eruptions!).  

All the files (student and instructor versions) of notebooks and homework questions are here:

https://github.com/gleesonm1/PetThermoTools/tree/master/docs/teaching_materials/VolatileSaturation
