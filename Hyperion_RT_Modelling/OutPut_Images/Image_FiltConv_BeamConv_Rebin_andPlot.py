from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from hyperion.util.constants import pc, c
from astropy.io import fits
import pyphot

from Image_FiltConv import Image_FilterConvolve
from Image_BeamConv import Image_BeamConvolve
from Image_Rebin import Image_Rebin


"""

Main image convolution script. 

Step 1: Image_FilterConvolve function - Convolving the four wavelength images output by RT modelling with the Filter profiles of their instruments.

Step 2:  Image_BeamConvolve function - Convolving the fileter convolved images with the beams of each instrument. 

"""

Source_Distance = 268 * pc #Units = pc
OutPutFiles = ["Model_FourShells_EvenMass.rtout", "Model_FourShells_UnEvenMass.rtout", "Model_InnerShell_Only.rtout", "Model_OuterShell_Only_K2010.rtout", "Model_OuterShell_Only_M2010.rtout", "Model_ShellThree_Only.rtout", "Model_ShellTwo_Only.rtout"] #RT output files
Filter_Library = pyphot.get_library(fname="Ampere_FiterProfile_Library.hdf5") #Getting Instrument Filters for Filter Convolution of the SED

#Filter names of the image wavelengths as given in the Ampere Filter porifles library. 
#'HERSCHEL_PACS_BLUE', 'HERSCHEL_PACS_RED', 'JCMT_SCUBA2_450', 'JCMT_SCUBA2_850'


#Needed for Beam Conv
model_max_envelope_size = 80 #Values used in RT modelling
model_arcsec_size = model_max_envelope_size * 3 #Values used in RT modelling
model_pix_size = 400 #Values used in RT modelling -  (400 x 400 pixels)
model_pix_arcsec = model_arcsec_size / model_pix_size #Size of each pixel in model image in arcsecs


for filename in OutPutFiles:

    model = ModelOutput(filename)

    wavelength_group = np.array([1, 2, 3, 4]) #wavelength group in data cube to exract the image at the wanted wavelength. Groups (in the order written out in RT modelling code) =>> 0 - SED (ignore here as we need the images) // 1 - PACS 70; 2 - PACS 160; 3 - SCUBA2 450; 4 - SCUBA2 850

    for group_val in wavelength_group: 

		#print(group_val)

		#Applying Filer convolution to the Model image at each wavelength
        Filter_Name, Filter_ConvImage = Image_FilterConvolve(model, Source_Distance, Filter_Library, group_val)

        #Creating .fits files of the filter convolved image at each wavelength
        hdu = fits.PrimaryHDU(Filter_ConvImage)
        hdulist = fits.HDUList([hdu])
        hdulist[0].header['Filter'] = Filter_Name[0]
        hdulist[0].header['BUNIT'] = "Jy"
        hdulist.writeto(filename.strip('.rtout') + '_' + Filter_Name[0] + '_FilterConv.fits', overwrite=True)

        #Applying Beam convolution onto to the Filter convolved image
        BeamConv_Image = Image_BeamConvolve(Source_Distance, group_val, Filter_Name[0], Filter_ConvImage, model_pix_arcsec)

        #Creating .fits files of the filter and then beam convolved image at each wavelength
        hdu = fits.PrimaryHDU(BeamConv_Image)
        hdulist = fits.HDUList([hdu])
        hdulist[0].header['Filter'] = Filter_Name[0]
        hdulist[0].header['BUNIT'] = "Jy"
        hdulist.writeto(filename.strip('.rtout') + '_' + Filter_Name[0] + '_FiltandBeamConv.fits', overwrite=True)

        #Now we want to rebin the image and convert it to Jy/arcsec**2
        Rebinned_image = Image_Rebin(BeamConv_Image, Filter_Name[0], model_pix_arcsec)
        
        #Creating .fits files of the Rebinned Beam Convolved Images - Final Output required to form Radial profiles. Image units = Jy/arcsec**2
        hdu = fits.PrimaryHDU(Rebinned_image)
        hdulist = fits.HDUList([hdu])
        hdulist[0].header['Filter'] = Filter_Name[0]
        hdulist[0].header['BUNIT'] = "Jy/arcsec**2"
        hdulist.writeto(filename.strip('.rtout') + '_' + Filter_Name[0] + '_FinalRebinnedImage.fits', overwrite=True)

        #plt.imshow(Rebinned_image)
        #plt.show()


				

		






##########################





















