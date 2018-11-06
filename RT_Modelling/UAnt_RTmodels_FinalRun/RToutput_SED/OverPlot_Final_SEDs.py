import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
#import pylatex
import matplotlib.ticker as mtick
import warnings


"""
Overplot all Final filter convolved SEDs of the different models with the Observed SED
"""

font = {'family' : 'normal',
        'size'   : 18,
	'weight' : 'medium'}

#Heavily Edited for U Ant!! Don't use for anything else!!!!!!!!!!

plt.rc('font', **font)


#star = 'uant'
#distance = 268.097  #Units = pc (Mc.Donald et al., 2017 from Hippacos distances from van Loon 2007??)


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


#Plotting
fig = plt.figure(figsize=(12, 10)) #(width,height)
ax = fig.add_subplot(1, 1, 1)

a = np.argsort(obs_wavelength) #Sort the arrays in ascending wavelength order to plot correctly.


# Plot SED for each inclination
ax.errorbar(obs_wavelength[a], outer_M2010_flux[a], fmt='-', color='red', markersize=6, capsize=2, label='Model_Outer_M2010') #Outer Shell (M2010) + Cont Emission - Shell size from Maercker+2010
ax.errorbar(obs_wavelength[a], outer_K2010_flux[a], fmt='-', color='blue', markersize=6, capsize=2, label='Model_Outer_K2010') #Outer Shell (K2010) + Cont Emission - Shell size from Kerschbaum+2010
ax.errorbar(obs_wavelength[a], inner_flux[a], fmt='-', color='green', markersize=6, capsize=2, label='Model_Inner') # Model SED Inner Shell + Cont Emission
ax.errorbar(obs_wavelength[a], fourshells_even_flux[a], fmt='-', color='purple', markersize=6, capsize=2, label='Model_FourShells_EvenMass') #Model SED Four Shells + Cont Emission - Dust Mass Evenly spread between four shells
ax.errorbar(obs_wavelength[a], fourshells_uneven_flux[a], fmt='-', color='orange', markersize=6, capsize=2, label='Model_FourShells_UnevenMass') #Model SED Four Shells + Cont Emission - Dust Mass Unevenly spread between four shells
ax.errorbar(obs_wavelength[a], cont_only_flux[a], fmt='-', color='violet', markersize=6, capsize=2, label='Model_NoShells') #Model SED Four Shells + Cont Emission - Dust Mass Evenly spread between four shells
ax.errorbar(obs_wavelength[a], obs_flux[a], yerr=obs_flux_unc[a], fmt='o', color='black', markersize=6, capsize=2, label='Observed SED') #Observed SED


ax.legend(fontsize=11, loc=1)
ax.set_xscale('log') #setting x axis to log scale
ax.set_yscale('log') #setting y axis to log scale
ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel(r'$F_\lambda$ [Jy]')
ax.set_xlim(1, 1500.)

SED_name = 'OverPlot_SEDs.png'
fig.savefig(SED_name)
plt.show()






