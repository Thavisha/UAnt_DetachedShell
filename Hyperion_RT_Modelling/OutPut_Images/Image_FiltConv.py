from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from hyperion.util.constants import pc, c
from astropy.io import fits
import pyphot

'''

Convolving the four wavelength images output by RT modelling with the Filter profiles of their instruments. This is step 1 of the image conolution process. Step 2 is the Beam convolution. 

'''

def Image_FilterConvolve(model, Source_Distance, Filter_Library, group_val):

	image = model.get_image(group=int(group_val), inclination=0, distance=Source_Distance, units='Jy') #Full Image data cube holding all four wavelength images		

	#Getting Frequency and Wavelength of the images in the output cubes (1 cube for each wavelength)
	RT_Image_Wavelengths = image.wav #image.wav = list of wavelengths created based on the limits given during image creation at each wavelengths in RT modelling stage.Pre Convolution
	RT_Image_Fluxes = image.val  #Fluxes of the image derived during RT modelling. Pre Convolution 

	#print('Wavelength of Image = ', RT_Image_Wavelengths, 'micron')
	#print(RT_Image_Fluxes.shape[:-1]) 

	Filter_ConvImage = np.zeros(RT_Image_Fluxes.shape[:-1]) #Create a empty array in the same shape as the 1st two (everything up to the 1 before last element) axes (400 x 400 in our case) to store the filter convolved photometry. 

	#Extracting the correct filters from the filt.library
	#Filter for PACS 70 image
	if RT_Image_Wavelengths[0] < 100: 
		Filter_Name = ['HERSCHEL_PACS_BLUE']
	
	#Filter for PACS 160 image
	elif RT_Image_Wavelengths[0] < 300:
		Filter_Name = ['HERSCHEL_PACS_RED']

	#Filter for SCUBA-2 450 image
	elif RT_Image_Wavelengths[0] < 500:
		Filter_Name = ['JCMT_SCUBA2_450']

	#Filter for SCUBA-2 850 image
	elif RT_Image_Wavelengths[0] < 1000:
		Filter_Name = ['JCMT_SCUBA2_850']


	#Extracting the correct filters from the filt.library	
	Filter = Filter_Library.load_filters(Filter_Name, interp=True, lamb=RT_Image_Wavelengths*pyphot.unit['micron'])  #pyphot.unit['micron'] tells pyphot that the wavlengths are in micron units.

	#Loop over each pixel in the image (400 x 400) to carry out filt.conv along the z axes
	#First get the row 
	for row in range(RT_Image_Fluxes.shape[0]): 
		#Then get the col. and therefore the pixel in the chosen row. 
		for col in range(RT_Image_Fluxes.shape[1]): 

			#Convolving with Filter Profiles
			filter_info, FiltConv_Image_Fluxes = pyphot.extractPhotometry(RT_Image_Wavelengths, RT_Image_Fluxes[row, col], Filter, Fnu=True, absFlux=False, progress=False) #filter_info saves info about the filt. profiles. SED flux units = Jy. #progress=False - won't print progress on terminal. 

			Filter_ConvImage[row, col] = FiltConv_Image_Fluxes

	#plt.imshow(Filter_ConvImage)
	#plt.show()

    

	return Filter_Name, Filter_ConvImage


