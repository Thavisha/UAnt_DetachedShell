import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import warnings


"""
Calculating the Reduced Chisquared values of each model with respect to the observed SED

When running this give an output filename in the command line to save the Chi Squared results in the file
eg: python RedChiSquared_SED.py > RedChiSquaredValues_SED.csv #This prints the output to the filename RedChiSquaredValues_SED.csv

"""


#Loading Observed SED of U Ant. 
obs_sed = Table.read('Observed_SED.csv', format='ascii.csv') 
obs_wavelength = obs_sed['Wavelength_micron']
obs_flux = obs_sed['Obs_SED_Flux_Jy']
obs_flux_unc = obs_sed['Obs_SED_Flux_Unc_Jy']

#Loading Convolved Model SED Outer Shell (M2010) + Cont Emission - Shell size from Maercker+2010
outer_M2010_sed = Table.read('Model_OuterShell_ContEm_M2010_Filter_Convolved_SED.csv', format='ascii.csv') 
outer_M2010_flux = outer_M2010_sed['Filt.Convolved_SED_Flux_Jy']


#Loading Convolved Model SED Outer Shell (K2010) + Cont Emission - Shell size from Kerschbaum+2010
outer_K2010_sed = Table.read('Model_OuterShell_ContEm_K2010_Filter_Convolved_SED.csv', format='ascii.csv') 
outer_K2010_flux = outer_K2010_sed['Filt.Convolved_SED_Flux_Jy']


#Loading Convolved Model SED Inner Shell + Cont Emission
inner_sed = Table.read('Model_InnerShell_ContEm_Filter_Convolved_SED.csv', format='ascii.csv') 
inner_flux = inner_sed['Filt.Convolved_SED_Flux_Jy']


#Loading Convolved Model SED Four Shells + Cont Emission - Dust Mass Evenly spread between four shells
fourshells_even_sed = Table.read('Model_FourShells+ContEm_EvenMass_Filter_Convolved_SED.csv', format='ascii.csv') 
fourshells_even_flux = fourshells_even_sed['Filt.Convolved_SED_Flux_Jy']


#Loading Convolved Model SED Four Shells + Cont Emission - Dust Mass Unevenly spread between four shells
fourshells_uneven_sed = Table.read('Model_FourShells_ContEm_UnevenMass_Filter_Convolved_SED.csv', format='ascii.csv') 
fourshells_uneven_flux = fourshells_uneven_sed['Filt.Convolved_SED_Flux_Jy']


#Loading Convolved Model SED No Shells - Only Cont emission
cont_only_sed = Table.read('Model_ContEm_Only_Filter_Convolved_SED.csv', format='ascii.csv') 
cont_only_flux = cont_only_sed['Filt.Convolved_SED_Flux_Jy']


##### ChiSqaured Generation #################

#Val = obs_wavelength >= 1.2
#print('Val =', Val)

print("Model Type, Reduced Chisqaured Result") #Prints a heading in the file given when running python. 

#"Reduced" Chisqaured (Divide Chisquared by the number points in obs) for Model SED Outer Shell (M2010) + Cont Emission - Shell size from Maercker+2010
outer_M2010_RedChiSq = ( np.sum ( ( (obs_flux - outer_M2010_flux) / obs_flux_unc) **2 ) ) / np.shape(obs_flux)[0]
#print(np.shape(obs_flux))
print("outer_M2010_RedChiSq,", outer_M2010_RedChiSq) 


#Chisqaured for Model SED Outer Shell (K2010) + Cont Emission - Shell size from Kerschbaum+2010
outer_K2010_RedChiSq = np.sum ( ( (obs_flux - outer_K2010_flux) / obs_flux_unc) **2 ) / np.shape(obs_flux)[0]
print("outer_K2010_RedChiSq,", outer_K2010_RedChiSq) 


#Chisqaured for Model SED Inner Shell + Cont Emission
inner_RedChiSq = np.sum ( ( (obs_flux - inner_flux) / obs_flux_unc) **2 ) / np.shape(obs_flux)[0]
print("inner_RedChiSq,", inner_RedChiSq) 


#Chisqaured for Model SED Four Shells + Cont Emission - Dust Mass Evenly spread between four shells
fourshells_even_RedChiSq = np.sum ( ( (obs_flux - fourshells_even_flux) / obs_flux_unc) **2 ) / np.shape(obs_flux)[0]
print("fourshells_even_RedChiSq,", fourshells_even_RedChiSq) 


#Chisqaured for Model SED Four Shells + Cont Emission - Dust Mass Unevenly spread between four shells
fourshells_uneven_RedChiSq = np.sum ( ( (obs_flux - fourshells_uneven_flux) / obs_flux_unc) **2 ) / np.shape(obs_flux)[0]
print("fourshells_uneven_RedChiSq,", fourshells_uneven_RedChiSq) 


#Chisqaured for Model SED No Shells - Only Cont emission
cont_only_RedChiSq = np.sum ( ( (obs_flux - cont_only_flux) / obs_flux_unc) **2 ) / np.shape(obs_flux)[0]
print("cont_only_RedChiSq,", cont_only_RedChiSq) 



























