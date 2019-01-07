# SAUCES

Stellar Activity for Understanding and Characterizing Exoplanetary Systems

For applications of this code see:
Roettenbacher & Vida (2018; https://arxiv.org/pdf/1810.04762.pdf)

These codes provide a way to investigate stellar activity, and they are currently written in IDL.

Included (so far!) 

* SAUCES_Kepler_CBV_remove.pro:  IDL process for removing cotrending basis vectors from Kepler light curves (DR25 CBVs).

* SAUCES_period_find_FFT.pro:  IDL process for finding the most likely period of a light curve and for identifying if a light curve is likely showing evolving spots or something more stable (e.g., eclipsing binaries or pulsations).  

* SAUCES_Kepler_local_min.pro:  IDL process for smoothing the light curves and finding the local minima

* SAUCES_Kepler_min_histplots.pro:  IDL process for making histograms of the local minima

* SAUCES_Kepler_missing_data.pro:  IDL process for identifying gaps in the data

* SAUCES_Kepler_flares.script:  IDL script for putting of these preparation processes together

* SAUCES_Kepler_flares_compare.pro:  IDL process for comparing the location of starspots and the times of flares

* Kepler_CBV_list.txt:  Text file with the names of the Kepler DR25 CBVs.

* Kepler_quarter_suffixes.txt:  Text file with the suffixes of the Kepler quarters.  

If you use or reference these codes, please cite Roettenbacher, R. M. & Vida, K. 2018, ApJ, 868, 3 (https://doi.org/10.3847/1538-4357/aae77e).
