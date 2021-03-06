import numpy as np
from hyperion.model import AnalyticalYSOModel
from hyperion.util.constants import rsun, lsun, au, msun, year, c, sigma #sigma = notation for stefan-boltzman constant


"""

Model Set up for a model Envelope with Four detached shells from thermal pulse (shell 4 from optical) 

Total shell Mass devided unevenly between the four shells with 50% in shell 4, 5% in shell 3, 22.5% in shell 2, 22.5% in shell 1 

"""

def FourShells_UnEvenMass(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature, stellar_radius, total_shell_mass, uant_distance, expansion_velocity,  MaxCSE_outer_rad, Photon_Number, CPU_Number):

	# Initalize the model
	model = AnalyticalYSOModel() #We use the YSO model with modified parameters to match an AGB star with a detached shell

	# Set the stellar parameters
	model.star.spectrum = (nu, fnu)
	model.star.luminosity = stellar_luminoisity
	#model.star.mass = 2 * msun #Guesstimate 1.3 - 3 Msun #Not needed for RT modelling normally. Very difficult to estimate.
	model.star.radius = stellar_radius



	#Power-law spherically symmetric envelope - Only Detached Shell - Shell 1 
	envelope_shell1 = model.add_power_law_envelope()
	envelope_shell1.mass = total_shell_mass * 0.225              # 22.5% of total dust in shell 1
	envelope_shell1.rmin = (23.5 * uant_distance) * au        #GD+2001&2003          
	envelope_shell1.rmax = (26.5 * uant_distance) * au       #GD+2001&2003   
	envelope_shell1.r_0 = envelope_shell1.rmin
	envelope_shell1.power = -2                        # Radial power #Constant Outflow envelope
	envelope_shell1.dust =  dust_file_name 

	#Power-law spherically symmetric envelope - Only Detached Shell - Shell 2  
	envelope_shell2 = model.add_power_law_envelope()
	envelope_shell2.mass = total_shell_mass * 0.225              # 22.5% of total dust in shell 2
	envelope_shell2.rmin = (34 * uant_distance) * au        #GD+2001&2003          
	envelope_shell2.rmax = (40 * uant_distance) * au        #GD+2001&2003  
	envelope_shell2.r_0 = envelope_shell2.rmin
	envelope_shell2.power = -2                        # Radial power #Constant Outflow envelope
	envelope_shell2.dust =  dust_file_name 

	#Power-law spherically symmetric envelope - Only Detached Shell - Shell 3 
	envelope_shell3 = model.add_power_law_envelope()
	envelope_shell3.mass = total_shell_mass * 0.05        # 5% of total dust in shell 3 - Gas dominated
	envelope_shell3.rmin = (42 * uant_distance) * au      #Mearcker+2010            
	envelope_shell3.rmax = (44 * uant_distance) * au     #Mearcker+2010     
	envelope_shell3.r_0 = envelope_shell3.rmin
	envelope_shell3.power = -2                        # Radial power #Constant Outflow envelope
	envelope_shell3.dust =  dust_file_name #Using a standard Hyperion Dust model. Needs to be modfied a lot for U Ant

	#Power-law spherically symmetric envelope - Only Detached Shell - Shell 4 
	envelope_shell4 = model.add_power_law_envelope()
	envelope_shell4.mass = total_shell_mass * 0.5              # 50% of total dust in shell 4
	envelope_shell4.rmin = (46.5 * uant_distance) * au          # Inner radius - Maercker+2010
	envelope_shell4.rmax = (53.5 * uant_distance) * au          # Outer radius - Maercker+2010
	envelope_shell4.r_0 = envelope_shell4.rmin
	envelope_shell4.power = -2                        # Radial power #Constant Outflow envelope
	envelope_shell4.dust =  dust_file_name 


	envelope_CSE_Max = (MaxCSE_outer_rad * uant_distance) * au    # Outer radius- estimating the total emission will be ~ this radius - pre&post thermal pulse. Units = cm


	# Set up grid to run the model on
	model.set_spherical_polar_grid_auto(100, 1, 1) #(n_r, n_theta, n_phi) - in a spherical envelope theta and phi are uniform and 1. n_r = number of r radii to seperate the grid to. 


	# Use raytracing to improve s/n of thermal/source emission
	model.set_raytracing(raytracing=True)

	# Set up SED - Get SED output
	sed = model.add_peeled_images(sed=True, image=False)
	sed.set_uncertainties(uncertainties=True)
	sed.set_viewing_angles([45], [45]) #(np.linspace(0., 90., 10), np.repeat(45., 10)) #Veiw the source at 45degrees from the pole and theta = 45
	sed.set_wavelength_range(100, 0.3, 1200.) #(n_wav, wav_min, wav_max) #Get the SED from 1micron to 2000micron in 100 wavelengths. 

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
	Input_FourShells_UnevenMass = 'Model_FourShells_UnEvenMass.rtin' #File which saves all the input above so that fortran can run internally for RT modelling
	Output_FourShells_UnevenMass = 'Model_FourShells_UnEvenMass.rtout'
	
	model.write(Input_FourShells_UnevenMass)
	model.run(Output_FourShells_UnevenMass, mpi=True, n_processes=CPU_Number)

	return Input_FourShells_UnevenMass, Output_FourShells_UnevenMass

















