from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from hyperion.util.constants import pc, c
from astropy.io import fits
import pyphot
from astropy.io import fits
from scipy import ndimage
from scipy import interpolate
from astropy.convolution import convolve, CustomKernel

"""
Convolving the fileter convolved images with the beams of each instrument. This is step 2 of image convolution. Step 1 was filter convolution

Bruteforce Convolution method based on https://www.pyimagesearch.com/2016/07/25/convolutions-with-opencv-and-python/ by Adrian Rosebrock. Last accessed 04/11/2018

"""



#Filter names of the image wavelengths as given in the Ampere Filter porifles library. 
#'HERSCHEL_PACS_BLUE', 'HERSCHEL_PACS_RED', 'JCMT_SCUBA2_450', 'JCMT_SCUBA2_850'

#PSF psf profile for SCUBA2
def Gaussian(x, fwhm, center=[0]):  #'b' in this case is the initial x axis before interpolating it all into one common x axis scale

          if len(x) == 1:
                    x = np.arange(0, size, 1, float)
   
          if center is None:
                    x0 = size // 2
          else:
                    x0 = center[0]
   
          psf_shape =  np.exp(-4*np.log(2) * ((x-x0) / fwhm)**2)  
          
          return psf_shape  #eqn from Dempsey et al.,2013 (SCUBA2 PSF gaussian psf derived from observations of Uranus)


#Bruteforce convolution
def Convolve(Filter_ConvImage, PSF_array):

    #Getting spatial dimensions of image and kernel
    (image_height, image_width) =  Filter_ConvImage.shape
    (kernel_height, kernel_width) =  PSF_array.shape

    #Pad border of image to ensure the kernel does not go outisde the image. 
    padding = (kernel_height - 1) // 2
    padded_Image = np.pad(Filter_ConvImage, ((padding, padding), (padding, padding)), mode="edge")
    #Padded_Image = opencv.copyMakeBorder(Filter_ConvImage, padding, padding, padding, padding, opencv.BORDER_REPLICATE)

    #Define empty array to hold output image for post conv. 
    BeamCov_Image = np.zeros((image_height, image_width))

    #Convolve each pixel of Image with the kernel
    for y in np.arange(padding, image_height + padding):
        for x in np.arange(padding, image_width + padding):

            #Extracting the region we want to convolve over (aroud the chosen pixel) and obtain it's center (i.e: the chosen pixel value)
            Conv_region = padded_Image[y - padding:y + padding + 1, x - padding:x + padding + 1]

            #Carrying out convolution by taking the element wise multiplication between Conv_region and the kernel followed by matrix summation (multiplies the two arrays then sums them together)
            conv_product = (Conv_region * PSF_array).sum()

            #Updating the empty output array with the convolved information
            BeamCov_Image[y - padding, x - padding] = conv_product

    return BeamCov_Image



#Carrying out Beam Conv, following Filter conv. Filter_Name, Filter_ConvImage = Output from Filter conv. 
def Image_BeamConvolve(Source_Distance, group_val, Filter_Name, Filter_ConvImage, model_pix_arcsec): 

    #Calling on the Beams
    if Filter_Name == 'HERSCHEL_PACS_BLUE':
        Beam_file = fits.open('PACS_PSF_70.fits')
        Beam_image = Beam_file[0].data[1000-200:1000+201,1000-200:1000+201]
        print(Beam_image.shape)
        Beam_image_pix_size = 1 #1pix = 1" in PACS 70 PSF image
        Beam_image_coord = np.indices(Beam_image.shape)    
        Beam_x_coord = Beam_image_coord[0] * Beam_image_pix_size #* Beam_image_pix_size - Converting the image coord from arcsec to pix size 
        Beam_y_coord = Beam_image_coord[1] * Beam_image_pix_size
        Interp_function = interpolate.interp2d(np.unique(Beam_x_coord), np.unique(Beam_y_coord), Beam_image, kind='cubic', fill_value=0., bounds_error=False)

        if np.int(Beam_image.shape[0]*(Beam_image_pix_size/model_pix_arcsec)) % 2 == 0:
            Upsampled_array = np.zeros((np.int(Beam_image.shape[0]*(Beam_image_pix_size/model_pix_arcsec))+1 , 1+np.int(Beam_image.shape[1]*(Beam_image_pix_size/model_pix_arcsec))))
        else:
            Upsampled_array = np.zeros((np.int(Beam_image.shape[0]*(Beam_image_pix_size/model_pix_arcsec)), np.int(Beam_image.shape[1]*(Beam_image_pix_size/model_pix_arcsec))))  #Creating an empty array to up sample the PSF to match the pix size of the RT modelled image. Size of the empty array is  larger than the original PSF fits file.
        Upsampled_array_coord =  np.indices(Upsampled_array.shape)    
        Upsampled_x_coord = Upsampled_array_coord[0] * model_pix_arcsec
        Upsampled_y_coord = Upsampled_array_coord[1] * model_pix_arcsec
        PSF_array = Interp_function(np.unique(Upsampled_x_coord), np.unique(Upsampled_y_coord))

    elif Filter_Name == 'HERSCHEL_PACS_RED':
        Beam_file = fits.open('PACS_PSF_160.fits')
        Beam_image = Beam_file[0].data[1000-200:1000+201,1000-200:1000+201]
        print(Beam_image.shape)
        Beam_image_pix_size = 1 #1pix = 1" in PACS 70 PSF image
        Beam_image_coord = np.indices(Beam_image.shape)    
        Beam_x_coord = Beam_image_coord[0] * Beam_image_pix_size #* Beam_image_pix_size - Converting the image coord from arcsec to pix size 
        Beam_y_coord = Beam_image_coord[1] * Beam_image_pix_size
        Interp_function = interpolate.interp2d(np.unique(Beam_x_coord), np.unique(Beam_y_coord), Beam_image, kind='cubic', fill_value=0., bounds_error=False)

        if np.int(Beam_image.shape[0]*(Beam_image_pix_size/model_pix_arcsec)) % 2 == 0:
            Upsampled_array = np.zeros((np.int(Beam_image.shape[0]*(Beam_image_pix_size/model_pix_arcsec))+1 , 1+np.int(Beam_image.shape[1]*(Beam_image_pix_size/model_pix_arcsec))))
        else:
            Upsampled_array = np.zeros((np.int(Beam_image.shape[0]*(Beam_image_pix_size/model_pix_arcsec)), np.int(Beam_image.shape[1]*(Beam_image_pix_size/model_pix_arcsec))))  #Creating an empty array to up sample the PSF to match the pix size of the RT modelled image. Size of the empty array is  larger than the original PSF fits file.
        Upsampled_array_coord =  np.indices(Upsampled_array.shape)    
        Upsampled_x_coord = Upsampled_array_coord[0] * model_pix_arcsec
        Upsampled_y_coord = Upsampled_array_coord[1] * model_pix_arcsec
        PSF_array = Interp_function(np.unique(Upsampled_x_coord), np.unique(Upsampled_y_coord))


    elif Filter_Name == 'JCMT_SCUBA2_450':
        fwhm = 7.9
        fwhm_sec = 25
        npix = int( (fwhm_sec*4) / model_pix_arcsec) #Number pixels in the psf kernel. Kernal size must be odd so that the center is placed in the center.
        if npix % 2 == 0: #%2 - find the remainder of dividing by 2. =0 check if remainder is 0 and therefore even or odd
            npix += 1 #If remainder is zero add 1
        PSF_array = np.zeros((npix , npix))
        indices = np.indices(PSF_array.shape) - npix//2    #256
        radii = np.sqrt(indices[0]**2 + indices[1]**2) * model_pix_arcsec
        PSF_array = (Gaussian(radii, fwhm=fwhm)*0.6) + (Gaussian(radii, fwhm=fwhm_sec)*0.4) #SCUBA-2 beam shape from Dempsey et al., 2013

    elif Filter_Name == 'JCMT_SCUBA2_850':
        fwhm = 13
        fwhm_sec = 48
        npix = int( (fwhm_sec*4) / model_pix_arcsec) #Number pixels in the psf kernel. Kernal size must be odd so that the center is placed in the center.
        if npix % 2 == 0: #%2 - find the remainder of dividing by 2. =0 check if remainder is 0 and therefore even or odd
            npix += 1 #If remainder is zero add 1
        PSF_array = np.zeros((npix , npix))
        indices = np.indices(PSF_array.shape) - npix//2    #256
        radii = np.sqrt(indices[0]**2 + indices[1]**2) *model_pix_arcsec
        PSF_array = (Gaussian(radii, fwhm=fwhm)*0.75) + (Gaussian(radii, fwhm=fwhm_sec)*0.25) #SCUBA-2 beam shape from Dempsey et al., 2013    

    print(PSF_array.shape,Filter_ConvImage.shape)
    psf_sum = np.sum(PSF_array)
    print(psf_sum)
    #BeamConv_Image = ndimage.convolve(Filter_ConvImage, PSF_array, mode='constant', cval=0.0) #Scipy Convolution - Didn't work!!
    #BeamConv_Image = convolve(Filter_ConvImage, CustomKernel(PSF_array)) #Astropy Convolution - Didn't work!!

    BeamConv_Image = Convolve(Filter_ConvImage, PSF_array/psf_sum) #Bruteforce Convolution. PSF_array/psf_sum = To normalise the PSF correctly

    return BeamConv_Image



    







































