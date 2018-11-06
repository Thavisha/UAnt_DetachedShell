from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
from scipy import interpolate
from scipy.integrate import simps
from numpy import trapz
import matplotlib.ticker as mtick
import warnings
from glob import glob


"""
Create Radial profiles from the final Rebinned Image of the RT models

Modified Rad prof code used for Obs. 

"""

#Radial Profile
def radial_profile(data, center):

    x, y = np.indices((data.shape))
    rad = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    rad = rad.astype(np.int)

    radialbin = np.bincount(rad.ravel(), data.ravel())
    n_rad = np.bincount(rad.ravel())
    calc_radprofile = radialbin / n_rad
    
    return calc_radprofile


#Defining input fits files and their parameters
def profile(filename):

    #Loading fits file
    fitsFile = fits.open(filename) 

    #PACS 70 micron
    if filename[-28:-24] == 'BLUE': #elif=else if

            img = fitsFile[0].data # Herschel PACS 70 RT model
            pix_scale = 1.6 # PACS 70 micron pixel scale -- 1pix=1.6"

            #Defining CENTRAL pixel of obj -- (Y, X) -- Give Y coord first
            cx =  img.shape[0] // 2 #central pixel x value 
            cy =  img.shape[1] // 2 #central pixel y value 


    #PACS 160 micron      	         	
    elif filename[-27:-24] == 'RED':      	

            img = fitsFile[0].data ## Herschel PACS 70 RT model
            pix_scale = 3.2 # PACS 160 micron pixel scale -- 1pix=3.2"

            #Defining CENTRAL pixel of obj -- (Y, X) -- Give Y coord first
            cx =  img.shape[0] // 2 #central pixel x value 
            cy =  img.shape[1] // 2 #central pixel y value 


    #SCUBA2 450 micron
    elif filename[-27:-24] == '450': # SCUBA2 450 RT model
        
        img = fitsFile[0].data 
        pix_scale = 2 # SCUBA2 450 micron pixel scale -- 1pix=2"

        #Defining CENTRAL pixel of obj -- (Y, X) -- Give Y coord first
        cx =  img.shape[0] // 2 #central pixel x value 
        cy =  img.shape[1] // 2 #central pixel y value 


    #SCUBA2 850 micron
    elif filename[-27:-24] == '850':# SCUBA2 850 RT model
        
        img = fitsFile[0].data
        pix_scale = 4 # SCUBA2 850 micron pixel scale -- 1pix=4"

        #Defining CENTRAL pixel of obj -- (Y, X) -- Give Y coord first
        cx =  img.shape[0] // 2 #central pixel x value 
        cy =  img.shape[1] // 2 #central pixel y value 


    else:
        raise NotImplementedError('No fits file found')


    center = (cy-1, cx-1)
    radprof = radial_profile(img, center)
    x = np.zeros_like(radprof)  
    for i in range(len(x)): 
        x[i]=i*pix_scale  #converting x axis from pixels to arcsec using the given pixel scales ('pix_scale')

    return radprof, x, center 




############################################################################################

################################# Generating the profiles ##################################

if __name__=="__main__":

    
    filenames = glob('*_FinalRebinnedImage.fits')
    print(filenames)
    #filenames=['Model_OuterShell_ContEm_K2010_HERSCHEL_PACS_BLUE_FinalRebinnedImage.fits','Model_OuterShell_ContEm_K2010_HERSCHEL_PACS_RED_FinalRebinnedImage.fits', 'Model_OuterShell_ContEm_K2010_JCMT_SCUBA2_450_FinalRebinnedImage.fits', 'Model_OuterShell_ContEm_K2010_JCMT_SCUBA2_850_FinalRebinnedImage.fits'] 


    for filename in filenames:

        radprof_obj, x_axis, center = profile(filename)

        warnings.filterwarnings("ignore") #Ignore any warnings


########################## Saving All Radial Profile data for each source so it can be used for plotting the radial and residual plots ############

        f = open(filename.strip('_FinalRebinnedImage.fits') + '_RadialProfile.csv', 'w')
        f.write("x_axis(arcsec)" + "," + "Stellar_Profile(Jy/arcsec^2)" + "\n")
        for j in range(len(x_axis)):
            outstring = str(x_axis[j]) + "," + str(radprof_obj[j]) + "\n"
            f.write(outstring)
        f.close()      	


