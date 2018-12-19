###################################################################################################################################
    ########### Radial Profile Generation ################################
###################################################################################################################################

- Script to plot the image of the observation at a chosen SCUBA-2 wavelength or HARP wavelength.

- The PSF shape and three sigma brightness extension is overplotted on the image.

- Script for both SCUBA-2 wavelegnths and HARP.  

** This script can be used for both 450 micron obs and 850 micron obs by changing the wavelength varaible at the top of the script. We have only generated the 850micron obs images for this paper. **


###### Input files required ####################
1) .csv table file containing source central X and Y positions and three sigma radius - Source_Information_SCUBA2_850.csv
2) .fits SCUBA-2/HARP observation of chosen source. 


######## Output ####################

- Image of the source with overplotted PSF & three sigma brightness extension. As well as the four shells observed in the optical. Axes provided as RA and Dec in degrees as well as the surface brightness colour scale in mJy/arcsec^2. 































