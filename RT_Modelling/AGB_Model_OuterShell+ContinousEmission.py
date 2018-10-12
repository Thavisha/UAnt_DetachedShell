import numpy as np
from hyperion.model import AnalyticalYSOModel
from hyperion.util.constants import rsun, lsun, au, msun, year, c, sigma #sigma = notation for stefan-boltzman constant



############## Model Set up for a model Envelope with Outermost detached shell from thermal pulse (shell 4 from optical) + continous Outflow with to account for pre&post thermal pulse emission ############################

#Parameters
uant_distance = 268 #pc
stellar_spectrum_name = 'Aringer+2009_StellarSpectrum.ascii' #File Cols = wav-Angstroms, fcont-NotNeeded-Ignore, nufnu-Convert to fnu ('kt04000g+3.5z-2.0.ascii')
dust_file_name = 'hyperion-dust-0.1.0/dust_files/kmh94_3.1_full.hdf5'
total_shell_mass = 2e-5 * msun #From my SED fit (point-point) - Fitted after CO subtraction from 850
present_day_mlr = 1.34e-10  #Present day MLR(DPR) dervid from GRAMS for the the NESS target selection # units = Msun/yr
expansion_velocity = (4.5 * 10e5) * year #CO expsansion vel averaged from Bands6&3 from Kerschbaum+2017 converted from 4.5 km/s to cm/year for easy conversion to derive total constant outflow mass later on. 



# Initalize the model
m = AnalyticalYSOModel() #We use the YSO model with modified parameters to match an AGB star with a detached shell

#Stellar Spectrum - YSO model for testing
#wav, fnu = np.loadtxt(stellar_spectrum_name, unpack=True) #YSO star - Needs to be replaced with a good C star model
#nu = c / (wav * 1.e-4)


# Read in stellar spectrum for Aringer+2009_StellarSpectrum.ascii ##### Does this need to be measured at the distance of U Ant like in Sundar's file ?????????????????????????????!!!!!!!!!!!!!!!!!!!
wav, fcont, nufnu = np.loadtxt(stellar_spectrum_name, unpack=True) #wav, fcont, nufnu. Hyperion needs nu and fnu so we convert out file values to this
nu = c / (wav * 1.e-8) #(wav * 1.e-8) to convert from angstrom to cm
fnu = nufnu / nu

extraploation_nu = 3e8 #Extrapoplate up to 3e8 Hz (=100 cm) to have the stellar spectrum at long wavelengths 
extrapolation_function = ((extraploation_nu/nu[-1]) ** 2 ) * fnu[-1]
fnu = np.append(fnu, extrapolation_function)
nu = np.append(nu, extraploation_nu)
print(nu, fnu)



# Set the stellar parameters #Taken from current version of U Ant paper using GRAMS output
m.star.spectrum = (nu, fnu)
m.star.luminosity = 6000 * lsun
#m.star.mass = 2 * msun #Guesstimate 1.3 - 3 Msun #Not needed for RT modelling normally. Very difficult to estimate.
m.star.radius = np.sqrt (m.star.luminosity / (4 * sigma * np.pi * (2800 ** 4) ) ) #To get stellar radius(R) - L=4PiR^2T^4 (T=2800 K) #radius output in cm


#Power-law spherically symmetric envelope - Only Detached Shell - Need to add constant outflow envelope on top of this to account for the ML after the thermal pulse. 
envelope_shell = m.add_power_law_envelope()
envelope_shell.mass = total_shell_mass          
envelope_shell.rmin = (46.5 * uant_distance) * au          # Inner radius converted to au - from Maercker+2010
envelope_shell.rmax = (53.5 * uant_distance) * au          # Outer radius converted to au - from Maercker+2010
envelope_shell.r_0 = envelope_shell.rmin
envelope_shell.power = -2                        # Radial power #Constant Outflow envelope
envelope_shell.dust =  dust_file_name #Using a standard Hyperion Dust model. Needs to be modfied a lot for U Ant

#Power-law spherically symmetric envelope - Constant outflow envelope on top of Shell envelope to account for the ML before&after the thermal pulse. 
envelope_current = m.add_power_law_envelope()
envelope_current.rmin = m.star.radius * 5                       # Inner radius - 5r_star converted to au - estimation #Units = cm
envelope_current.rmax = (70 * uant_distance) * au          # Outer radius 70" - estimating the total emission will be ~70" - pre&post thermal pulse. Units = cm
envelope_current.r_0 = envelope_current.rmin
envelope_current.power = -2  # Radial power #Constant Outflow envelope
envelope_current.mass = ( ( (envelope_current.rmax - envelope_current.rmin) / expansion_velocity ) * present_day_mlr ) * msun #Derive total constant outflow mass from present day mlr and expansion vel
envelope_current.dust = dust_file_name  #Using a standard Hyperion Dust model. Needs to be modfied a lot for U Ant



# Set up grid
m.set_spherical_polar_grid_auto(100, 1, 1) #(n_r, n_theta, n_phi) - in a spherical shell theta and phi are uniform and 1



# Use raytracing to improve s/n of thermal/source emission
m.set_raytracing(raytracing=True)

# Set up SED - Get SED output
sed = m.add_peeled_images(sed=True, image=False)
sed.set_uncertainties(uncertainties=True)
sed.set_viewing_angles([45], [45]) #(np.linspace(0., 90., 10), np.repeat(45., 10)) #Veiw the source at 45degrees from the pole and theta = 45
sed.set_wavelength_range(100, 0.3, 1200.) #(n_wav, wav_min, wav_max) #Get the SED from 1micron to 2000micron in 100 wavelengths

#Set up Image - Get image output - RT calculated for bins of wavelengths
image = m.add_peeled_images(sed=False, image=True)
image.set_uncertainties(uncertainties=True)
image.set_viewing_angles([45], [45]) #(np.linspace(0., 90., 10), np.repeat(45., 10)) #Veiw the source at 45degrees from the pole and theta = 45
image.set_wavelength_range(1, 845, 855) #(n_wav, wav_min, wav_max) - Only get the image at between 800-900 micron for Radial profiling (scuba-2 filter profile range 750-1000micron) 
image.set_image_size(400, 400) #Size of image in pixels
image.set_image_limits(-1.5 * envelope_current.rmax, 1.5 * envelope_current.rmax, -1.5 * envelope_current.rmax, 1.5 * envelope_current.rmax)


#Set up Monochromatic Images - Get Monochromatic image output - RT calculated for each wavelength in given array
#image = m.add_peeled_images(sed=False, image=True)
#wavelength_array = np.arange(800,910, 10) #get an array from 800 to 900 in 10 steps
#image.set_monochromatic(monochromatic=True, wavelengths=wavelength_array)
#image.set_uncertainties(uncertainties=True)
#image.set_viewing_angles([45], [45]) #(np.linspace(0., 90., 10), np.repeat(45., 10)) #Veiw the source at 45degrees from the pole and theta = 45
#image.set_wavelength_range(1, 845, 855) #(n_wav, wav_min, wav_max) - Only get the image at between 800-900 micron for Radial profiling (scuba-2 filter profile range 750-1000micron) 
#image.set_image_size(400, 400) #Size of image in pixels
#image.set_image_limits(-1.5 * envelope_current.rmax, 1.5 * envelope_current.rmax, -1.5 * envelope_current.rmax, 1.5 * envelope_current.rmax)







# Set number of photons
m.set_n_photons(initial=1e7, imaging=1e7, raytracing_sources=1e7, raytracing_dust=1e7)

# Set number of temperature iterations and convergence criterion
m.set_n_initial_iterations(5)
m.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)

# Write out file
m.write('AGB_Model_Example_5.rtin')
m.run('AGB_Model_Example_5.rtout', mpi=True)

















