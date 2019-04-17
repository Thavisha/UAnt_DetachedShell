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

#Model 1: Inner shell only  (All shell dust mass put into shell 1 defined by G2001&2003)
inner_sed = Table.read('Model_InnerShell_Only_Filter_Convolved_SED.csv', format='ascii.csv') 
inner_flux = inner_sed['Filt.Convolved_SED_Flux_Jy']

#Model 2: Shell Two Only (All shell dust mass put into shell 2 defined by G2001&2003)
shellTwo_sed = Table.read('Model_ShellTwo_Only_Filter_Convolved_SED.csv', format='ascii.csv') 
shellTwo_flux = shellTwo_sed['Filt.Convolved_SED_Flux_Jy']

#Model 3: Shell Three Only (All shell dust mass put into shell 3 defined by G2001&2003 and M2010)
shellThree_sed = Table.read('Model_ShellThree_Only_Filter_Convolved_SED.csv', format='ascii.csv') 
shellThree_flux = shellThree_sed['Filt.Convolved_SED_Flux_Jy']

#Model 4: Shell 4 Only with Radius of shell as given in M2010 - All dust mass in shell 4
outer_M2010_sed = Table.read('Model_OuterShell_Only_M2010_Filter_Convolved_SED.csv', format='ascii.csv') 
outer_M2010_flux = outer_M2010_sed['Filt.Convolved_SED_Flux_Jy']

#Model 5: Shell 4 Only with Radius of shell as given in K2010 - All dust mass in shell 4
outer_K2010_sed = Table.read('Model_OuterShell_Only_K2010_Filter_Convolved_SED.csv', format='ascii.csv') 
outer_K2010_flux = outer_K2010_sed['Filt.Convolved_SED_Flux_Jy']

#Model 6: Four shells Model - Shell mass distributed unevenly based on literature values
fourshells_uneven_sed = Table.read('Model_FourShells_UnEvenMass_Filter_Convolved_SED.csv', format='ascii.csv') 
fourshells_uneven_flux = fourshells_uneven_sed['Filt.Convolved_SED_Flux_Jy']

#Model 7: Four Shells Model - Shell Mass distributed Evenly
fourshells_even_sed = Table.read('Model_FourShells_EvenMass_Filter_Convolved_SED.csv', format='ascii.csv') 
fourshells_even_flux = fourshells_even_sed['Filt.Convolved_SED_Flux_Jy']




##### REDUCED (Divide Chisquared by the number points in obs) ChiSqaured Generation #################

#Val = obs_wavelength >= 1.2
#print('Val =', Val)

print("Model Type, Reduced Chisqaured Result") #Prints a heading in the file given when running python. 

#Reduced Chisqaured for Model SED Inner Shell Only
inner_RedChiSq = np.sum ( ( (obs_flux - inner_flux) / obs_flux_unc) **2 ) / np.shape(obs_flux)[0]
print("Model_InnerShell_Only,", inner_RedChiSq) 

#Reduced Chisqaured for Model SED Shell Two Only
shellTwo_RedChiSq = np.sum ( ( (obs_flux - shellTwo_flux) / obs_flux_unc) **2 ) / np.shape(obs_flux)[0]
print("Model_ShellTwo_Only,", shellTwo_RedChiSq) 

#Reduced Chisqaured for Model SED Shell Two Only
shellThree_RedChiSq = np.sum ( ( (obs_flux - shellThree_flux) / obs_flux_unc) **2 ) / np.shape(obs_flux)[0]
print("Model_ShellThree_Only,", shellThree_RedChiSq) 

#Reduced Chisqaured for Model SED Outer Shell (M2010) Shell size from Maercker+2010
outer_M2010_RedChiSq = ( np.sum ( ( (obs_flux - outer_M2010_flux) / obs_flux_unc) **2 ) ) / np.shape(obs_flux)[0]
#print(np.shape(obs_flux))
print("Model_OuterShell_Only_M2010,", outer_M2010_RedChiSq) 


#Reduced Chisqaured for Model SED Outer Shell (K2010)  Shell size from Kerschbaum+2010
outer_K2010_RedChiSq = ( np.sum ( ( (obs_flux - outer_K2010_flux) / obs_flux_unc) **2 ) ) / np.shape(obs_flux)[0]
#print(np.shape(obs_flux))
print("Model_OuterShell_Only_K2010,", outer_K2010_RedChiSq) 


#Reduced  Chisqaured for Model SED Four Shells - Dust Mass Unevenly spread between four shells
fourshells_uneven_RedChiSq = np.sum ( ( (obs_flux - fourshells_uneven_flux) / obs_flux_unc) **2 ) / np.shape(obs_flux)[0]
print("Model_FourShells_UnEvenMass,", fourshells_uneven_RedChiSq) 

#Reduced  Chisqaured for Model SED Four Shells - Dust Mass Evenly spread between four shells
fourshells_even_RedChiSq = np.sum ( ( (obs_flux - fourshells_even_flux) / obs_flux_unc) **2 ) / np.shape(obs_flux)[0]
print("Model_FourShells_EvenMass,", fourshells_even_RedChiSq) 































