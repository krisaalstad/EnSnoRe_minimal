# Lakes_minimal

This repository hopefully contains the minimum amount of data and MATLAB code that will enable you to reproduce the snow reanalysis over the Mammoth Lakes basin that were conducted in Aalstad et al. (2020).

To do so, run these scripts in the following order:

1) Move to the Forcing directory and run gen_TS_forcing.m to downscale the NLDAS forcing in "all_NLDAS.mat".
2) Run agg_forcing.m to aggregate the forcing to the daily timestep of the model.
3) Move to the Simulation directory and run main.m to run the snow reanalysis for 2016-2019 (takes roughly 3-4 hours). 

Feel free to then compare the results generated in Simulation/result_20XX.mat to the independent ASO validation data stored in the Validation/Gridded_Val directory.

References: Aalstad, K. et al. (2020): Snow reanalysis using multiple optical satellite constellations, About to be submitted. 

N.B. I have only included the processed land-cover, DEM, satellite, and ASO data. The raw data is not included since it is too large to store on github. Furthermore, to keep the code base manageable I have not included some of the pre-processing scripts for the land-cover, DEM, satellite, AVIRIS, and ASO data. However, if you are interested in these steps feel free to get in touch with me. 

Cheers,
Kris

