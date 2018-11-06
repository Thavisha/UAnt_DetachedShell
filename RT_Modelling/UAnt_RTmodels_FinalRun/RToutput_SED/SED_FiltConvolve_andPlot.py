import matplotlib.pyplot as plt
import numpy as np
from hyperion.model import ModelOutput
from hyperion.util.constants import pc
import pyphot

Source_Distance = 268 * pc
OutPutFiles = ['Model_OuterShell_ContEm_M2010.rtout', 'Model_ContEm_Only.rtout', 'Model_InnerShell_ContEm.rtout', 'Model_FourShells+ContEm_EvenMass.rtout', 'Model_FourShells_ContEm_UnevenMass.rtout', 'Model_OuterShell_ContEm_K2010.rtout'] #RT output files

#Getting Instrument Filters for Filter Convolution of the SED
Filter_Library = pyphot.get_library(fname="Ampere_FiterProfile_Library.hdf5")

#Filter names of the photometric points in my SED as given in the Ampere Filter porifles library. !!!!!!!!! DIRBES STUFF ADD !!!!!!!!!!
Filter_Names = np.array(['2MASS_J', '2MASS_H', '2MASS_Ks', 'COBE_DIRBE_1.25', 'COBE_DIRBE_2.2', 'COBE_DIRBE_3.5', 'COBE_DIRBE_4.9', 'AKARI_S9W', 'AKARI_L18W', 'AKARI_FIS_N60', 'AKARI_FIS_WIDES', 'AKARI_FIS_WIDEL', 'IRAS_12', 'IRAS_25', 'IRAS_60', 'IRAS_100', 'HERSCHEL_PACS_BLUE', 'HERSCHEL_PACS_RED', 'HERSCHEL_SPIRE_PSW', 'HERSCHEL_SPIRE_PMW', 'HERSCHEL_SPIRE_PLW', 'JCMT_SCUBA2_450', 'JCMT_SCUBA2_850'])


for filename in OutPutFiles:

	model = ModelOutput(filename)
	sed = model.get_sed(group=0, inclination='all', aperture=-1, distance=Source_Distance, units='Jy')

	#print(sed.wav)
	#print(sed.val)
	
	RT_SED_wavelengths = sed.wav #sed.wav = list of wavelengths created based on the limits (100, 0.3, 1200.) given during SED creation in RT modelling stage. Pre Convolution with Filters
	RT_SED_Fluxes = sed.val[0] #Fluxes of the SED derived during RT modelling. Pre Convolution with Filters #[0] - Hyperion has gives a list of arrays instead of an array. So we need to pick the zeroth element array even if the rest are empty.

	#Extracting the correct filters from the filt.library	
	Filters = Filter_Library.load_filters(Filter_Names, interp=True, lamb=RT_SED_wavelengths*pyphot.unit['micron'])  #pyphot.unit['micron'] tells pyphot that the wavlengths are in micron units. 

	#Convolving with Filter Profiles
	filter_info, Convolved_Model_SED_Fluxes = pyphot.extractPhotometry(RT_SED_wavelengths, RT_SED_Fluxes, Filters, Fnu=True, absFlux=False)#, progress=False) #filter_info saves info about the filt. profiles. SED flux units = Jy

	Conv_SED_Wavelengths = np.array([a.magnitude for a in filter_info]) #Extracting the convolved SED wavelengths from the filter_info. Units=micron
	
	
	#Final Convolved SED Model Output Table File
	Table_name = filename.strip('.rtout') + '_Filter_Convolved_SED.csv' #filename.strip('.rtout') - Get rid of .rtout from filename and add the +.. string after for the table name
	f = open(Table_name, 'w')
	f.write("Filter_Name" + "," + "Wavelength_micron" + "," + "Filt.Convolved_SED_Flux_Jy" + "\n")
	f.close()

	f = open(Table_name, 'a')
	for i in range(len(Filter_Names)):
		outstring = str(Filter_Names[i]) + ","  + str(Conv_SED_Wavelengths[i]) + ","  + str(Convolved_Model_SED_Fluxes[i]) + "\n"
		f.write(outstring)
	f.close()


	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)

	# Plot SED for each inclination
	ax.errorbar(Conv_SED_Wavelengths, Convolved_Model_SED_Fluxes, fmt='o', color='black') #yerr=sed.unc - Can you convolve error bars???

	ax.set_xscale('log') #setting x axis to log scale
	ax.set_yscale('log') #setting y axis to log scale
	ax.set_xlabel(r'$\lambda$ [$\mu$m]')
	ax.set_ylabel(r'$F_\nu$ [Jy]')
	ax.set_xlim(0.1, 1500.)

	SED_name = filename.strip('.rtout') + '_Filter_Convolved_SED.png'

	fig.savefig(SED_name)

	plt.show()
















