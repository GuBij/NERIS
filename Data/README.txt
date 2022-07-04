The simulated dose rates are provided in the folder Data/doseRates. The employed models are
-	Advection-diffusion particle model with forest parameterization (F)
-	Advection-diffusion particle model with open field parameterization (OF)
-	Langevin particle model with forest parameterization (LF)
-	Langevin particle model with open field parameterization (LOF)
-	Gaussian plume model with Bultynck-Malet parameterization (GBM)
-	Gaussian plume model with Pasquill-Gifford parameterization (GPG)
The suffixes '60m' and '10m' in the folder name represent the used release height.
The SI folder contains the simulated dose rates by the LF model used for Table 3.7 p.64.
The Monin-Obukhov lengts in the files MOST_Forest.txt and MOST_OF.txt in the Data folder have been generated using the executable genMeteo.

The dose rates with the open field (OF) and forest (F) parameterization in the folder Data/doseRates have been generated using the executable simDose.
It is shown in the folder Data/doseRates/inputFiles how the run directory should be set up for the dose rate simulation:
-	The proper set-up of the run directory is obtained by running 'NERIS_forest.sh dirList.forest' for the forest and 'NERIS_OF.sh dirList.OF' 
	for the open field parameterization
-	If the advection-diffusion model instead of the Langevin model needs to be run, then remove '$option' after the call to genMeteo in the just mentioned sh files
-	The parameter set-up for the NERIS experiment is provided in the file 'inputDict' in the '.orig' folders in Data/doseRates/inputFiles
-	The file 'script.pbs' in the '.orig' folders shows how the code can be run in parallel

The structure of the meteofiles 'meteo<day><month>.txt' in the folder Data/doseRates/inputFiles/meteo is as follows
	<temperature at 8 m (°C)>	<temperature at 114 m (°C)>		<wind speed (m/s)>		<wind direction (°)>	<elevation (°)>		<\sigma_azi (°)>	<\sigma_elev (°)>
with \sigma the standard devation of the horizontal (azi) or the vertical wind direction (elev).
