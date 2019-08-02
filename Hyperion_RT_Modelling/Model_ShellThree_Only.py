import numpy as np
from hyperion.model import AnalyticalYSOModel
from hyperion.util.constants import rsun, lsun, au, msun, year, c, sigma #sigma = notation for stefan-boltzman constant


"""

Model Set up for a model Envelope with 3rd detached shell (shell 3 from optical). All the shell mass is put in shell 3 in this case. 

"""

def ShellThree_Only(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature,  stellar_radius, total_shell_mass, uant_distance, expansion_velocity, MaxCSE_outer_rad, Photon_Number, CPU_Number):

	# Initalize the model
	model = AnalyticalYSOModel() #We use the YSO model with modified parameters to match an AGB star with a detached shell

	# Set the stellar parameters
	model.star.spectrum = (nu, fnu)
	model.star.luminosity = stellar_luminoisity
	#model.star.mass = 2 * msun #Guesstimate 1.3 - 3 Msun #Not needed for RT modelling normally. Very difficult to estimate.
	model.star.radius = stellar_radius


	#Power-law spherically symmetric envelope - Only Detached Shell - Need to add constant outflow envelope on top of this to account for the ML after the thermal pulse. 
	envelope_shell = model.add_power_law_envelope()
	envelope_shell.mass = total_shell_mass          
	envelope_shell.rmin = (41.5 * uant_distance) * au          # Shell 3 inner radius converted to au - from GD+2001&2003 AND M2010
	envelope_shell.rmax = (44.5 * uant_distance) * au          # Shell 3 outer radius converted to au - from GD+2001&2003 AND M2010
	envelope_shell.r_0 = envelope_shell.rmin
	envelope_shell.power = -2                        # Radial power #Constant Outflow envelope
	envelope_shell.dust =  dust_file_name #Using a standard Hyperion Dust model. Needs to be modified a lot for U Ant

	envelope_CSE_Max = (MaxCSE_outer_rad * uant_distance) * au    # Outer radius- estimating the total emission will be ~ this radius - pre&post thermal pulse. Units = cm


	# Set up grid to run the model on
	model.set_spherical_polar_grid_auto(100, 1, 1) #(n_r, n_theta, n_phi) - in a spherical envelope theta and phi are uniform and 1. n_r = number of r radii to seperate the grid to. 


	# Use raytracing to improve s/n of thermal/source emission
	model.set_raytracing(raytracing=True)

	# Set up SED - Get SED output
	sed = model.add_peeled_images(sed=True, image=False)
	sed.set_uncertainties(uncertainties=True)
	sed.set_viewing_angles([45], [45]) #(np.linspace(0., 90., 10), np.repeat(45., 10)) #Veiw the source at 45degrees from the pole and theta = 45
	sed.set_wavelength_range(100, 0.3, 1200.) #(n_wav, wav_min, wav_max) #Get the SED from 1micron to 2000micron in 100 wavelengths

	#Set up Image - Get image output - RT calculated for bins of wavelengths - PACS 70
	image_70 = model.add_peeled_images(sed=False, image=True)
	image_70.set_uncertainties(uncertainties=True)
	image_70.set_viewing_angles([45], [45]) #(np.linspace(0., 90., 10), np.repeat(45., 10)) #Veiw the source at 45degrees from the pole and theta = 45
	image_70.set_wavelength_range(30, 60, 90) #(n_wav, wav_min, wav_max) - PACS 70 filter profile limits 60 - 90 and get 30 image cube to get images at ~1 micron per image 
	image_70.set_image_size(400, 400) #Size of image in pixels
	image_70.set_image_limits(-1.5 * envelope_CSE_Max, 1.5 * envelope_CSE_Max, -1.5 * envelope_CSE_Max, 1.5 * envelope_CSE_Max)

	#Set up Image - Get image output - RT calculated for bins of wavelengths - PACS 160
	image_160 = model.add_peeled_images(sed=False, image=True)
	image_160.set_uncertainties(uncertainties=True)
	image_160.set_viewing_angles([45], [45]) #(np.linspace(0., 90., 10), np.repeat(45., 10)) #Veiw the source at 45degrees from the pole and theta = 45
	image_160.set_wavelength_range(30, 130, 220) #(n_wav, wav_min, wav_max) - PACS 160 filter profile limits 130 - 220 and get 30 image cube to get images at ~3 micron per image  
	image_160.set_image_size(400, 400) #Size of image in pixels
	image_160.set_image_limits(-1.5 * envelope_CSE_Max, 1.5 * envelope_CSE_Max, -1.5 * envelope_CSE_Max, 1.5 * envelope_CSE_Max)

	#Set up Image - Get image output - RT calculated for bins of wavelengths - SCUBA-2 450
	image_450 = model.add_peeled_images(sed=False, image=True)
	image_450.set_uncertainties(uncertainties=True)
	image_450.set_viewing_angles([45], [45]) #(np.linspace(0., 90., 10), np.repeat(45., 10)) #Veiw the source at 45degrees from the pole and theta = 45
	image_450.set_wavelength_range(30, 410, 480) #(n_wav, wav_min, wav_max) -  - SCUBA-2 450 filter profile limits 790 - 940 and get 30 image cube to get images at ~2 micron per image 
	image_450.set_image_size(400, 400) #Size of image in pixels
	image_450.set_image_limits(-1.5 * envelope_CSE_Max, 1.5 * envelope_CSE_Max, -1.5 * envelope_CSE_Max, 1.5 * envelope_CSE_Max)

	#Set up Image - Get image output - RT calculated for bins of wavelengths - SCUBA-2 850
	image_850 = model.add_peeled_images(sed=False, image=True)
	image_850.set_uncertainties(uncertainties=True)
	image_850.set_viewing_angles([45], [45]) #(np.linspace(0., 90., 10), np.repeat(45., 10)) #Veiw the source at 45degrees from the pole and theta = 45
	image_850.set_wavelength_range(30, 790, 940) #(n_wav, wav_min, wav_max) - SCUBA-2 850 filter profile limits 790 - 940 and get 30 image cube to get images at ~5 micron per image 
	image_850.set_image_size(400, 400) #Size of image in pixels
	image_850.set_image_limits(-1.5 * envelope_CSE_Max, 1.5 * envelope_CSE_Max, -1.5 * envelope_CSE_Max, 1.5 * envelope_CSE_Max)


	# Set number of photons
	model.set_n_photons(initial=Photon_Number, imaging=Photon_Number, raytracing_sources=Photon_Number, raytracing_dust=Photon_Number)

	# Set number of temperature iterations and convergence criterion
	model.set_n_initial_iterations(5)
	model.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)

	# Write out file
	Input_ShellThree_Only = 'Model_ShellThree_Only.rtin' #File which saves all the input above so that fortran can run internally for RT modelling
	Output_ShellThree_Only = 'Model_ShellThree_Only.rtout'
	
	model.write(Input_ShellThree_Only)
	model.run(Output_ShellThree_Only, mpi=True, n_processes=CPU_Number)

	return Input_ShellThree_Only, Output_ShellThree_Only

















