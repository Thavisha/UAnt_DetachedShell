from astropy.io import fits
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from hyperion.util.constants import pc, c
from astropy.io import fits



Source_Distance = 268 #Units = pc
RT_Output_Filename = 'AGB_Model_Example_6.rtout'
Fits_File_Name = 'AGB_Model_Example_6.fits'
#FitsCube_File_Name = 'AGB_Model_Example_Cube_6.fits'

print(c)



# Retrieve image cube as before
m = ModelOutput(RT_Output_Filename)
image = m.get_image(group=1, inclination=0, distance=Source_Distance * pc, units='Jy') #group 1 instead of 0. #0 is where the SED info is held. #Distance units-pc therefore we need to multiply by the pc to convert the pc from pc to cm??????????????????????
plt.imshow(image.val[:,:,0])
#plt.show()

#Getting Frequency and Wavelength of the output image
Freq = image.nu
Wavelength = (c*(10**4)) / Freq #Units=micron. #c is originally in cm/s so we multiply by 10^4 to convert it to micron/s to get wav in micron

print('Frequency of Image = ', Freq, 'Hz')
print('Wavelength of Image = ', Wavelength, 'micron')


#Plotting Fits file
# The image extracted above is a 3D array. We can write it out to FITS. # We need to swap some of the directions around so as to be able to use. # the ds9 slider to change the wavelength of the image.
#fits.writeto(FitsCube_File_Name, image.val.swapaxes(0, 2).swapaxes(1, 2), overwrite=True) #Only if your plotting an array of wavelengths and need a cube. No need for a single wavelength image

# We can also just output one of the wavelengths
fits.writeto(Fits_File_Name, image.val[:, :, 0], overwrite=True)


#Adding Wavelength of the fits file to the header. 
hdulist = fits.open(Fits_File_Name)
hdulist[0].header["Wavelen"] = float(Wavelength), 'micron'
hdulist[0].header["Filter"] = 'corresponding model for JCMT/SCUBA-2 850 micron'
hdulist.writeto('AGB_Model_Example_5.fits', overwrite=True)

































