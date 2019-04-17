############################### Radiative Transfer Modelling using Hyperion ############################################

- Carry out Hyperion RT modelling in order to test multiple shell scenarios for U Ant which are then compared to the observed Radial and global SED profiles. 

- Qualitatively compare the re-sulting surface brightness profiles to those derived from theSCUBA-2 and PACS observation. The best fits model SEDsand surface brightness profiles allow us to narrow down the most likely CSE structure of U Ant. Quantitative comparisons made via The chi-squaredvalues  per  observed  data  point  (χ^2_p)

- Hyperion - Robitaille T. P., 2011, A&A, 536, A79

- The  input  stellar  parameters  are derived by fitting  the  observed SED of U Ant from optical to the mid-IR to the COMARCS model grid (Aringer et al., 2009) 

- The script is written with each of the tested model seperated into a single python funtion file and then all called together with the parent python script (Run_UAnt_Models.py). Therefore you require the function files for all six tested models and the parent script to run the full set of models. The parent script can be edited to run only a single or few models depending on the user's needs. 

- The output SEDs and Synthetic observations are covolved with the respective filter profiles and beams. 


###### Input files required ####################

1) Ampere_FiterProfile_Library.hdf5 - This file is too large to be uploaded to Github. It will be made public with the AMPERE introudction paper by Scicluna et al., in prep. In the mean time please email the lead-author for the filter library.  
2) Cstar_DustProperties.hdf5 Once again this file is too large to be uploaded to Github. Please email the lead-author for this. The scripts and required input files required to build this is provided in the repository.  
3) Aringer+2009_StellarSpectrum.ascii

