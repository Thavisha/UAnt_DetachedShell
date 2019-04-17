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

c_cm = 29979245800.0 #Speed of light in cm/s

#Loading Observed SED of U Ant. 
obs_sed = Table.read('Observed_SED.csv', format='ascii.csv') 
obs_wavelength = obs_sed['Wavelength_micron']
obs_flux = obs_sed['Obs_SED_Flux_Jy']
obs_flux_unc = obs_sed['Obs_SED_Flux_Unc_Jy']


#Best Fit stellar photosphere model from COMARCS grid 
# Read in stellar spectrum for Aringer+2009_StellarSpectrum.ascii
wav, fcont, nuLnu = np.loadtxt("Aringer+2009_StellarSpectrum.ascii", unpack=True) #wav, fcont, nuLnu. 
nu = c_cm / (wav * 1.e-8) #(wav * 1.e-8) to convert from angstrom to cm
Lnu = nuLnu / nu #In reality Aringer Model gives nuLnu but the difference between that and nuFnu is scaling factor which Hyperion automatically account for by scaling up to the given stellar luminosity so we don't need worry about converting from one to the other

extraploation_nu = 3e11 #Extrapoplate up to 3e8 Hz (3e8=100 cm) (3e11=1mm)to have the stellar spectrum at long wavelengths 
extrapolation_function = ((extraploation_nu/nu[-1]) ** 2 ) * Lnu[-1]
Lnu = np.append(Lnu, extrapolation_function)
nu = np.append(nu, extraploation_nu)
wav_extrap_micron = (c_cm / nu) * (1e4) #* (10^4) to convert wavelength from cm to micron
#print(nu, Lnu)

#Scaling Factor for Stellar spectrum to scale it down to the observed SED (This done automatically in Hyperion for RT modelling but for plotting we need to do it again)
Band_2MassK_wav = 2.2024757594239204 #Wavelength of 2MASS_Ks band
Stellar_Spectrum_2MASSK_Flux_index = (np.abs(wav_extrap_micron - Band_2MassK_wav)).argmin()
Stellar_Spectrum_2MASSK_Flux = Lnu[Stellar_Spectrum_2MASSK_Flux_index]
Obs_SED_Flux_Kband_index = (np.abs(obs_wavelength - Band_2MassK_wav)).argmin()
Obs_SED_Flux_Kband = obs_flux[Obs_SED_Flux_Kband_index]
Scalling_Factor = Stellar_Spectrum_2MASSK_Flux / Obs_SED_Flux_Kband



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







#Plotting
fig = plt.figure(figsize=(12, 10)) #(width,height)
ax = fig.add_subplot(1, 1, 1)

a = np.argsort(obs_wavelength) #Sort the arrays in ascending wavelength order to plot correctly.


# Plot SED for each inclination

#Stellar Photosphere
ax.errorbar(wav_extrap_micron, (Lnu/Scalling_Factor), fmt='--', color='dimgrey', alpha=0.8, label='COMARCS model') #Stellar Photosphere

#Models
ax.errorbar(obs_wavelength[a], inner_flux[a], fmt='-', color='green', markersize=6, linewidth=1, capsize=2, label='Mss1') # Model SED Inner Shell Only
ax.errorbar(obs_wavelength[a], shellTwo_flux[a], fmt='-', color='violet', markersize=6, linewidth=1, capsize=2, label='Mss2') # Model SED Shell Two Only
ax.errorbar(obs_wavelength[a], shellThree_flux[a], fmt='-', color='springgreen', markersize=6, linewidth=1, capsize=2, label='Mss3') # Model SED Shell Three Only
ax.errorbar(obs_wavelength[a], outer_M2010_flux[a], fmt='-', color='red', markersize=6, capsize=2, linewidth=1, label='Mss4-M2010') #Outer Shell (M2010) Shell size from Maercker+2010
ax.errorbar(obs_wavelength[a], outer_K2010_flux[a], fmt='-', color='gold', markersize=6, capsize=2, linewidth=1, label='Mss4-K2010') #Outer Shell (K2010) Shell size from Kerschbaum+2010
ax.errorbar(obs_wavelength[a], fourshells_uneven_flux[a], fmt='-', color='blue', markersize=6, capsize=2, linewidth=1, label='Mfourshells') #Model SED Four Shells Dust Mass UnEvenly spread between four shells
#ax.errorbar(obs_wavelength[a], fourshells_even_flux[a], fmt='-', color='purple', markersize=6, capsize=2, linewidth=1, label='Model_FourShells_EvenMass') #Model SED Four Shells - Dust Mass Evenly spread between four shells


#Observed SED
ax.errorbar(obs_wavelength[a], obs_flux[a], yerr=obs_flux_unc[a], fmt='o', color='black', markersize=6, capsize=2, label='Observed SED') #Observed SED

#Plotting SCUBA-2 450 data points in different Colour
b = obs_wavelength == 450
ax.errorbar(obs_wavelength[b], obs_flux[b], yerr=obs_flux_unc[b], fmt='s', color='black', markerfacecolor='lime', markersize=8, capsize=6, markeredgecolor='black', label='SCUBA-2 points')

#Legend put here to not have two scuba-2 labels.
ax.legend(fontsize=13.5, loc=1)

b = obs_wavelength == 850
ax.errorbar(obs_wavelength[b], obs_flux[b], yerr=obs_flux_unc[b], fmt='s', color='black', markerfacecolor='lime', markersize=8, capsize=6, markeredgecolor='black')



ax.set_xscale('log') #setting x axis to log scale
ax.set_yscale('log') #setting y axis to log scale
ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel(r'$F_\nu$ [Jy]')
ax.set_xlim(0, 1500.)

SED_name = 'OverPlot_SEDs.png'
fig.savefig(SED_name, bbox_inches='tight')
plt.show()






