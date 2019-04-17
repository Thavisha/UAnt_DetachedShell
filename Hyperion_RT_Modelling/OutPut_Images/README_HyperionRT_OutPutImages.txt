############################### Radiative Transfer Modelling using Hyperion - Output Images ############################################

Convolve the output Images of each model from running Hyperion with the filter profiles and the beam profiles at the four observed wavelengths.

Derive radial profiles from the filter and beam convolved synthetic observations.  

Carry out comparison between the observed Radial profiles and the model radial profiles using Reduced Chi sqaured per data point.

Image_FiltConv_BeamConv_Rebin_andPlot.py runs the filter and beam convolution as the parent script combining the seperate function scripts provided. 


###### Input files required ####################

1) .rtout files generated from Hyperion output for each of the six (or user chosen) models - The .rtout files generated for this paper are once again too large to be uploaded for GitHub. Please email lead author if they are required. 
2) Ampere_FiterProfile_Library.hdf5 - This file is too large to be uploaded to Github. It will be made public with the AMPERE introudction paper by Scicluna et al., in prep. In the mean time please email the lead-author for the filter library.  
3) Observed radial (.csv files) profiles derived in this paper (e.g: uant_UnInterpolated_850.csv)

