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

#Heavily Edited for U Ant!! Don't use for anything else!!!!!!!!!!


#Radial Profile
def radial_profile(data, center, uncertainty=False):
    #print np.indices((data.shape))
    x, y = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    if uncertainty==True: 
	radialprofile = tbin / (nr**(3./2)) #the decimal place after the 3 indicate the resulting value of 3/2 is a decimal value
	radialprofile = np.sqrt(radialprofile)
    else:
	radialprofile = tbin / nr
    return radialprofile
    
    
#Defining input fits files and their parameters
def profile(filename,pacs_70_x,pacs_70_y,pacs_160_x,pacs_160_y,scuba2_450_x,scuba2_450_y,scuba2_850_x,scuba2_850_y,uncertainty=False):

	#Loading fits file
    	fitsFile = fits.open(filename) 
    	
   	
    	if uncertainty==True: #uncertainty frames
        	try:
            		img = fitsFile['stDev'].data**2 #PACS error frame. **2=To convert into the uncertainty (combined with line 21) --radialprofile = tbin / (nr**(3./2))--
        	except:
            		try:
                		img = fitsFile['VARIANCE'].data[0] #SCUAB2 error frame. 
            		except:
                		raise NotImplementedError('Incorrect name for error file') #gives warning that file name for error file is incorrect
        
        
        #PACS 70 micron        		
	if filename[-7:-5] == '70': #elif=else if
		if filename[-11:-5] == 'PSF_70': #PACS PSF 70 micron -- PACS 70 PSF/Beam generated using Asteroid Vesta + Mars measurements - Bocchio et al 2016	
			img = fitsFile[0].data # Herschel PACS PSF data
			
			if uncertainty==True: #uncertainty frames
        			try:
            				img = fitsFile['stDev'].data**2 #PACS error frame. **2=To convert into the uncertainty (combined with line 21) --radialprofile = tbin / (nr**(3./2))--
        			except:
            				try:
                				img = fitsFile['VARIANCE'].data[0] #SCUAB2 error frame. 
            				except:
                				raise NotImplementedError('Incorrect name for error file') #gives warning that file name for error file is incorrect
				
			pix_scale = 1 # PACS PSF 70 micron pixel scale -- 1pix=1"
			img = (img/(pix_scale**2)) # to converting FLUX (y axis,entire image) from Jy/Pix to Jy/acrsec**2 
			
			#Defining CENTRAL pixel of PACS PSF -- (Y, X) -- Give Y coord first
		   	cx = 1000 #central pixel x value
		   	cy = 1000 #central pixel y value

		else:
		   	img = fitsFile['image'].data # Herschel PACS object 70 micron data
		   	
		   	if uncertainty==True: #uncertainty frames
        			try:
            				img = fitsFile['stDev'].data**2 #PACS error frame. **2=To convert into the uncertainty (combined with line 21) --radialprofile = tbin / (nr**(3./2))--
        			except:
            				try:
                				img = fitsFile['VARIANCE'].data[0] #SCUAB2 error frame. 
            				except:
                				raise NotImplementedError('Incorrect name for error file') #gives warning that file name for error file is incorrect
		   	
		   	pix_scale = 1.6 # PACS 70 micron pixel scale -- 1pix=1.6"
		   	img = (img/(pix_scale**2)) # to converting FLUX (y axis,entire image) from Jy/Pix to Jy/acrsec**2.
		   	
		   	#Defining CENTRAL pixel of obj -- (Y, X) -- Give Y coord first
		   	cx = pacs_70_x #central pixel x value #!!!CHANGE TO MATCH OBJECT!!!!
		   	cy = pacs_70_y #central pixel y value #!!!CHANGE TO MATCH OBJECT!!!!
		 		   	
			   	
			
	#PACS 160 micron	   	
	elif filename[-8:-5] == '160':	
			   	
		if filename[-12:-5] == 'PSF_160': #PACS PSF 160 micron -- PACS 160 PSF/Beam generated using Asteroid Vesta + Mars measurements - Bocchio et al 2016	
			img = fitsFile[0].data # Herschel PACS PSF data
			
			if uncertainty==True: #uncertainty frames
        			try:
            				img = fitsFile['stDev'].data**2 #PACS error frame. **2=To convert into the uncertainty (combined with line 21) --radialprofile = tbin / (nr**(3./2))--
        			except:
            				try:
                				img = fitsFile['VARIANCE'].data[0] #SCUAB2 error frame. 
            				except:
                				raise NotImplementedError('Incorrect name for error file') #gives warning that file name for error file is incorrect
			
			pix_scale = 1 # PACS PSF 160 micron pixel scale -- 1pix=1"
			img = (img/(pix_scale**2)) # to converting FLUX (y axis,entire image) from Jy/Pix to Jy/acrsec**2   
			
			#Defining CENTRAL pixel of obj -- (Y, X) -- Give Y coord first
		   	cx = 1000 #central pixel x value
		   	cy = 1000 #central pixel y value
			
		else:	
			img = fitsFile['image'].data # Herschel PACS object 160 micron data
			
			if uncertainty==True: #uncertainty frames
        			try:
            				img = fitsFile['stDev'].data**2 #PACS error frame. **2=To convert into the uncertainty (combined with line 21) --radialprofile = tbin / (nr**(3./2))--
        			except:
            				try:
                				img = fitsFile['VARIANCE'].data[0] #SCUAB2 error frame. 
            				except:
                				raise NotImplementedError('Incorrect name for error file') #gives warning that file name for error file is incorrect
			
			pix_scale = 3.2 # PACS 160 micron pixel scale -- 1pix=3.2"
		   	img = (img/(pix_scale**2)) # to converting FLUX (y axis,entire image) from Jy/Pix to Jy/acrsec**2.
		   	
		   	#Defining CENTRAL pixel of obj -- (Y, X) -- Give Y coord first
		   	cx = pacs_160_x  #central pixel x value #!!!CHANGE TO MATCH OBJECT!!!!
		   	cy = pacs_160_y #central pixel y value #!!!CHANGE TO MATCH OBJECT!!!!
		   	
		   	
	#SCUBA2 450 micron
	elif filename[-8:-5] == '450': # JCMT SCUBA2 object 450 micron data
			img = fitsFile[0].data[0] #use 0th frame (primary image/frame) of the fits file for SCUBA2
			
			if uncertainty==True: #uncertainty frames
        			try:
            				img = fitsFile['stDev'].data**2 #PACS error frame. **2=To convert into the uncertainty (combined with line 21) --radialprofile = tbin / (nr**(3./2))--
        			except:
            				try:
                				img = fitsFile['VARIANCE'].data[0] #SCUAB2 error frame. 
            				except:
                				raise NotImplementedError('Incorrect name for error file') #gives warning that file name for error file is incorrect
			
           		img[np.isnan(img)] = 0
           		
           		pix_scale = 2 # SCUBA2 450 micron pixel scale -- 1pix=2"
           		
           		#Defining CENTRAL pixel of obj -- (Y, X) -- Give Y coord first
		   	cx = scuba2_450_x #central pixel x value #!!!CHANGE TO MATCH OBJECT!!!!
		   	cy = scuba2_450_y #central pixel y value #!!!CHANGE TO MATCH OBJECT!!!!
		   	
           		
	#SCUBA2 850 micron
	elif filename[-8:-5] == '850': #JCMT SCUBA2 object 850 micron data
			img = fitsFile[0].data[0] #use 0th frame (primary image/frame) of the fits file for SCUBA2
			
			if uncertainty==True: #uncertainty frames
        			try:
            				img = fitsFile['stDev'].data**2 #PACS error frame. **2=To convert into the uncertainty (combined with line 21) --radialprofile = tbin / (nr**(3./2))--
        			except:
            				try:
                				img = fitsFile['VARIANCE'].data[0] #SCUAB2 error frame. 
            				except:
                				raise NotImplementedError('Incorrect name for error file') #gives warning that file name for error file is incorrect
			
           		img[np.isnan(img)] = 0
           		
           		pix_scale = 4 # SCUBA2 850 micron pixel scale -- 1pix=4"
           		
           		#Defining CENTRAL pixel of obj -- (Y, X) -- Give Y coord first
		   	cx = scuba2_850_x  #central pixel x value #!!!CHANGE TO MATCH OBJECT!!!!
		   	cy = scuba2_850_y  #central pixel y value #!!!CHANGE TO MATCH OBJECT!!!!
		   	
		   		
		
	    			
	else:
		raise NotImplementedError('No fits file found')
	    			
	
	center = (cy-1, cx-1)
	radprof = radial_profile(img, center, uncertainty)
	x = np.zeros_like(radprof)  
	for i in range(len(x)): 
 		x[i]=i*pix_scale  #converting x axis from pixels to arcsec using the given pixel scales ('pix_scale')
 		
	return radprof, x, center #????? ADD sigma here for PACS 160?????????



#PSF psf profile for SCUBA2
def Gaussian(x, fwhm, center=[0]):  #'b' in this case is the initial x axis before interpolating it all into one common x axis scale

     	'''
     	Calculates a gaussian function given a FWHM and x values

	Inputs:
     	x -
     	fwhm -
     	center -

    	Outputs:
    
	'''

	if len(x) == 1:
		x = np.arange(0, size, 1, float)
   
	if center is None:
		x0 = size // 2
	else:
		x0 = center[0]
   
	psf_shape =  np.exp(-4*np.log(2) * ((x-x0) / fwhm)**2)  
	
	return psf_shape, x  #eqn from Dempsey et al.2013 (SCUBA2 PSF gaussian psf derived from observations of Uranus)
		



############################################################################################

################################# Generating the profiles ##################################

if __name__=="__main__":



	star = 'uant' #!!!CHANGE TO MATCH OBJECT!!!!

	pacs_70_x = 495
	pacs_70_y = 506

	pacs_160_x = 250
	pacs_160_y = 253

	scuba2_450_x = 268    
	scuba2_450_y = 243
 
	scuba2_850_x = 131
	scuba2_850_y = 122
 

	filenames=[star+'_70.fits', star+'_160.fits', star+'_450.fits', star+'_850.fits'] #fits file names !!!CHANGE TO MATCH OBJECT!!!!
	#filenames=['epaqr_70.fits', 'epaqr_160.fits', 'epaqr_450.fits', 'epaqr_850.fits'] #fits file names !!!CHANGE TO MATCH OBJECT!!!!
	herschel_psf_files = ['PACS_PSF_70.fits', 'PACS_PSF_160.fits']
	
	x_interp = np.zeros(36) #36*4 = 144 approximate max x axis length needed (to interpolate onto the SCUBA2 850 pixel scale which is the lowest resolution we have) .
	for i in range(len(x_interp)): 
 		x_interp[i]=i*4 #setting up x axis for all wavelenghts to match the SCUBA2 850 micron data (to interpolate onto the SCUBA2 850 pixel scale which is the lowest resolution we have) 
 	
 	
 		
 	#Building the arrays to hold the data
 	stellar_profile = np.zeros((36,4)) #generate a blank/only zeros array to store the object profiles. #to loop over all 14 objct = np.zeros((36,4,14))
 	psf_profile = np.zeros((36,4)) #generate a blank/only zeros array to store the psf profiles.
 	unc_profile = np.zeros((36,4)) #generate a blank/only zeros array to store the unc profiles.
 	wavelengths = np.array([70,160,450,850]) #sq brackets = python list
 	fractional_unc_psf = np.zeros((4))
 	
	
	for filename in filenames: 
	
		
		if filename[-7:-5] == '70' or filename[-8:-5] == '160': #PACS fits files
                        

			#Object Radial Profile
                        radprof_obj,x_temp,center = profile(filename,pacs_70_x,pacs_70_y,pacs_160_x,pacs_160_y,scuba2_450_x,scuba2_450_y,scuba2_850_x,scuba2_850_y) #x_temp = temp x axis onto which individual profiles is built on to before the x axis's are all intrpltd to the same x axis so that the profiles can be comapred. 
                        radprof_obj_unc,x_temp,center = profile(filename,pacs_70_x,pacs_70_y,pacs_160_x,pacs_160_y,scuba2_450_x,scuba2_450_y,scuba2_850_x,scuba2_850_y,uncertainty=True) #error bars
                                              
                        if filename[-7:-5] == '70': 
                        	i_lam = 0 #defining the array col in which the profile will be filled in  to
                        else: 
                        	i_lam = 1 #defining the array col in which the profile will be filled in  to
			
			#PSF Profile
                        psffile = herschel_psf_files[i_lam] #i_lam in this case is the actual numbers 0 and 1 here for the 0th filename and 1st filename in herschel_psf_files
                        radprof_psf,x_psf,center = profile(psffile,pacs_70_x,pacs_70_y,pacs_160_x,pacs_160_y,scuba2_450_x,scuba2_450_y,scuba2_850_x,scuba2_850_y) #PSF(Beam) (0th frame). x_psf = equivalent of x_temp for radial profiles.
                       
                        
                        #Interpolating PSF profile. Each individual psf is interpolated onto their corresposnind radial profile temp x axis'
 			interp_psf = interpolate.interp1d(x_psf, radprof_psf, kind='cubic', fill_value=0., bounds_error=False) 
 			radprof_interp_psf = interp_psf(x_interp) #interpolated psf porifle 
			radprof_psf = interp_psf(x_temp) #psf's are interpolated onto the uninterpolated radial profile x axis. i.e: each individual psf is interpolated onto their corresposnind radial profile temp x axis's. 

		elif filename[-8:-5] == '450' or filename[-8:-5] == '850': #SCUBA2 fits files 
			
			#Obejct Radial Profile
			radprof_obj,x_temp,center = profile(filename,pacs_70_x,pacs_70_y,pacs_160_x,pacs_160_y,scuba2_450_x,scuba2_450_y,scuba2_850_x,scuba2_850_y) 
			radprof_obj_unc,x_temp,center = profile(filename,pacs_70_x,pacs_70_y,pacs_160_x,pacs_160_y,scuba2_450_x,scuba2_450_y,scuba2_850_x,scuba2_850_y,uncertainty=True) #error bars
			
			#PSF Profile. eqn from Dempsey et al. 2013
			if filename[-8:-5] == '450':
			
				alpha = 0.94 #value from Dempsey et al. 2013
				beta = 0.06  #value from Dempsey et al. 2013
				radprof_psf = (Gaussian(x_temp, fwhm=7.9)[0]*alpha) + (Gaussian(x_temp, fwhm=25)[0]*beta) #uniterploated scuba2 psf
				radprof_interp_psf = (Gaussian(x_interp, fwhm=7.9)[0]*alpha) + (Gaussian(x_interp, fwhm=25)[0]*beta) #interpolated SCUBA2 psfs
				i_lam = 2 #defining the array col in which the profile will be filled in  to
				
			else: 
				alpha = 0.98 #value from Dempsey et al. 2013
				beta = 0.02  #value from Dempsey et al. 2013
				radprof_psf = (Gaussian(x_temp, fwhm=13)[0]*alpha) + (Gaussian(x_temp, fwhm=48)[0]*beta) #uniterploated scuba2 psf
				radprof_interp_psf = (Gaussian(x_interp, fwhm=13)[0]*alpha) + (Gaussian(x_interp, fwhm=48)[0]*beta) #interpolated SCUBA2 psfs
				i_lam = 3 #defining the array col in which the profile will be filled in to
				
		


		#Scaling un-interpolated PSF porfile to match un-interpolated obejct profile (raising/lowering the uniterpolated PSF profile to match the object profile's peak y value (at x=zero)) 
		scimax=np.max(radprof_obj) #factor needed for creating the scale factor
		psf_scalefactor = scimax/np.max(radprof_psf)
		#print ('PSF Scalefactor = ', psf_scalefactor)
		radprof_psf_scaled = scimax/np.max(radprof_psf)*radprof_psf #scale factor = "scimax/np.max(radprof_psf)" = factor to which the psf is to be raised/lowered to match the interpolated radial profiles' peak y value at x=zero. 

	
		#Interpolating object. Interpolating all object radial profiles onto one common x axis (x_interp) so that that they can be compared side by side and used for SEDs
		interp_obj = interpolate.interp1d(x_temp, radprof_obj, kind='cubic', fill_value=0., bounds_error=False)
		radprof_interp_obj = interp_obj(x_interp)		
		scimax_interp=np.max(radprof_interp_obj) #factor needed for creating the scale factor.
                     
			
		#Interpolating unc's to match the interpolated object radial profiles
		interp_obj_unc = interpolate.interp1d(x_temp, radprof_obj_unc, kind='cubic', fill_value=0., bounds_error=False)
		radprof_interp_obj_unc = interp_obj_unc(x_interp)
		
		
		#Scaling interpolated Beam to match interpolated Object profile	
		psf_ineterp_scalefactor = scimax_interp/np.max(radprof_interp_psf)
		psf_interp_scaled = scimax_interp/np.max(radprof_interp_psf)*radprof_interp_psf #scale factor = "scimax_interp/np.max(radprof_interp_psf)" = factor to which the psf is to be raised/lowered to match the interpolated radial profiles' peak y value at x=zero. 
 		
 		
 		
			
 		
 		
		#Array building. Arrays to hold the profiles we generate. These arrays will then be used to be read into SED plotting.   
		#Add array col/row names
		stellar_profile[:,i_lam] = radprof_interp_obj
		psf_profile[:,i_lam] = psf_interp_scaled
		unc_profile[:,i_lam] =  radprof_interp_obj_unc
		fractional_unc_psf[i_lam] = radprof_obj_unc[0]/radprof_obj[0] #Fractional unc taken from peak (x=0) point of rad profile. For the PSF unc
		
		
		
		#Lists for plotting the UN-interpolated data. Here we're first trying to append a list called x.append in the try command. Then as that won't work since there is no x list created yet we're then creating an empty list x[] and then appending it to fill in a loop for each i. so the loop will go over for each i (i.e: wavelength)
		
		warnings.filterwarnings("ignore") #Ignore any warnings
		
		try: 
			x.append(x_temp)
		except: 
			x = []
			x.append(x_temp)
			
		
		try: 
			object_prof.append(radprof_obj)
		except:
			object_prof = []
			object_prof.append(radprof_obj)
			
		try: 
			psf_prof.append(radprof_psf_scaled)
		except:
			psf_prof = []
			psf_prof.append(radprof_psf_scaled)
			
		try: 
			object_prof_unc.append(radprof_obj_unc)
		except:
			object_prof_unc = []
			object_prof_unc.append(radprof_obj_unc)
			
			
			
############################### Secondary Profiles FOR INTERPOLATED PROFILES ############################	
		
	#Residual profile --> (Object - Beam)
	res_interp = stellar_profile - psf_profile
	res_interp_unc = np.sqrt((unc_profile**2)+((fractional_unc_psf*psf_profile)**2))


	#Fractional residual profile --> (Object - Beam)/Beam --> Residual as a fractional excess compared to Beam
	fracres_interp = (res_interp)/psf_profile #Fractional residual profile
	fracres_interp_unc = fracres_interp* np.sqrt(
						     ((res_interp_unc/res_interp)**2) +
						     (((fractional_unc_psf*psf_profile)/psf_profile)**2)
						     )


	#Printing data with the radial profile data for SED fitting in the terminal. #!!!THIS DOES WORK!!!! res_interp is printed as [70 data, 160 data, 450 data, 850 data]....
	#print filename 
	#print x_interp
	#print res_interp
	#print res_interp_unc
	


######## Calculating Results from the Radial Profiles for Interpolated Profiles ###################

	#Peak (Characteristic) Extension --> Determined from the from the Residual Profile
	for i in range(4):
	
		warnings.filterwarnings("ignore") #Ignore any warnings
		
		#Characteristic Extension (x value corresponing to peak y value)
		CharExt = x_interp[np.argmax(res_interp[:,i])] #print the x value which corresponds to the max y(=res) value
		print ('Characteristic Extension for Interpolated', wavelengths[i], 'um =', CharExt)  
	

		#Maximum Extension (Max radius at 3 sigma detection)
		SNR = res_interp[:,i]/res_interp_unc[:,i] #signal to noise ratio for 3sigma detection (3sigma detection=SNR)
		MaxExt = x_interp[SNR>=3] #print all x values SNR is greater than 3 (=3sigma) 
		#print ('Maximum Extension for Interpolated', wavelengths[i], 'um =', MaxExt[-1])
		try: print ('Maximum Extension for Interpolated', wavelengths[i], 'um =', MaxExt[-1]) 
		except:
			print ('Maximum Extension for Interpolated = No Extension') #when there is no extension for the source 
			#pass #sometimes there is no extension at all (mostly for SCUBA2 450 data) and then there is no value index -1 to be read out so an error will come up and the code will stop. to prevent this we add a try-except block where it tries to print the max ext and if it's not there it just moves on to the next wavelenght. 

		
		#x_interp = x_interp < float(MaxExt[-1])
		#z = [92, 68, 12, 44]
		#p = x_interp[i] < z[i]

	 	#Total Flux of the entire object (Flux under object radial profile)
		TotFlux = simps((2*np.pi*x_interp*stellar_profile[:,i]), x_interp)
		# print ('Total Flux for Interpolated', wavelengths[i], 'um =', TotFlux)

		#Flux under the PSF (flux under psf radial profile)
		PsfFlux = simps((2*np.pi*x_interp*psf_profile[:,i]), x_interp)
		#print ('PSF Flux for Interpolated', wavelengths[i], 'um =', PsfFlux)

		#Flux with extended portion (Total flux - PSF flux) 
		ExtFlux = TotFlux - PsfFlux 
		#print ('Extended Flux for Interpolated', wavelengths[i], 'um =', ExtFlux) 
	
		


############################# Secondary Profiles FOR UN-INTERPOLATED DATA ##################################### 

	
	for i in range(len(x)):
	
		warnings.filterwarnings("ignore") #Ignore any warnings
			
		#Residual profile --> (Object - Beam)
		try: 
			res.append(object_prof[i] - psf_prof[i])
			res_unc.append(np.sqrt((object_prof_unc[i]**2)+((fractional_unc_psf[i]*psf_prof[i])**2)))
		except: 
			res = []
			res_unc = []
			res.append(object_prof[i] - psf_prof[i])
			res_unc.append(np.sqrt((object_prof_unc[i]**2)+((fractional_unc_psf[i]*psf_prof[i])**2)))
			
			#Here we're first trying to append a list called res in the try command. Then as that won't work since there cont..
			#cont.. is no res list created we're then creating an empty list res[] and then appending it to fill in a cont..
			#cont.. loop for each i. so the loop will go over for each i (i.e: wavelength)
			
		#Fractional residual profile --> (Object - Beam)/Beam --> Residual as a fractional excess compared to Beam
		try:
			fracres.append(res[i]/psf_prof[i])
			fracres_unc.append(fracres[i]* np.sqrt(
						     ((res_unc[i]/res[i])**2) +
						     (((fractional_unc_psf[i]*psf_prof[i])/psf_prof[i])**2)
						     ))
			
		except:
			fracres = []
			fracres_unc = []
			fracres.append(res[i]/psf_prof[i])
			fracres_unc.append(fracres[i]* np.sqrt(
						     ((res_unc[i]/res[i])**2) +
						     (((fractional_unc_psf[i]*psf_prof[i])/psf_prof[i])**2)
						     ))
		
			
######## Calculating Results from the Radial Profiles for Un-iterpolated data ################### CHECK FLUX CODE!!?????????????????????

	#Peak (Characteristic) Extension --> Determined from the from the Residual Profile
	for i in range(4):
		
		warnings.filterwarnings("ignore") #Ignore any warnings
		
		
		#Peak (Characteristic) Length --> Determined from the from the Residual Profile
		CharExt = x[i][np.argmax(res[i])] #print the x value which corresponds to the max y(=rad_profile_res) value
		print ('Characteristic Extension for UN-Interpolated', wavelengths[i], 'um =', CharExt)
		

		#Max Extension (maximum amount of extended emission) --> Determined from the from the Residual Profile
		SNR = res[i]/res_unc[i] #signal to noise ratio for 3sigma detection (3sigma detection=SNR)
		MaxExt = x[i][SNR>=3] #print all x values SNR is greater than 3 (=3sigma)

		try: 
			#print ('Maximum Extension for UN-Interpolated', wavelengths[i], 'um =', MaxExt[-1])
			print ('Maximum Extension for UN-Interpolated', wavelengths[i], 'um =', MaxExt[0:])  
		
		except:
			print ('Maximum Extension for UN-Interpolated = NO Extension') #when there is no extension for the source 
			#pass #sometimes there is no extension at all (mostly for SCUBA2 450 data) and then there is no value index -1 to be read out so an error will come up and the code will stop. to prevent this we add a try-except block where it tries to print the max ext and if it's not there it just moves on to the next wavelenght. 
		
		#z = x[i] =< MaxExt[i]
		#b = x[i] < length[i]
		#x[i][b]
		#z = x[i] < str(MaxExt[-1])

		
		#Total Flux
		#TotFlux = simps((2*np.pi*x[z]*object_prof[z]), x[z])
		#print ('Total Flux for UN-Interpolated', wavelengths[i], 'um =', TotFlux)

		#Flux under the PSF 
		#PsfFlux = simps((2*np.pi*x[i][z]*psf_prof[i][z]),x[i][z])
		#print ('PSF Flux for UN-Interpolated', wavelengths[i], 'um =', PsfFlux)

		#Flux within extended portion 
		#ExtFlux = TotFlux - PsfFlux 
		#print ('Extended Flux for UN-Interpolated', wavelengths[i], 'um =', ExtFlux)
 
            
                
        
####################### Plotting the profile for all four wavelenghts from UN-interpolated data. #####################
#It's better plotting uniterpolated data since it looks better.


	for i in range(len(x)):
		fig = plt.figure(figsize=(10,11)) #(width,height)
		gs = gridspec.GridSpec(3, 1, height_ratios=[1.5, 1, 1]) #Note capital G and S in .GridSpec
		gs.update(hspace=0)
		#length = [100, 60, 30, 60] #x axis arcseconds lengths for 70, 160, 450, 850 in this specific order!
		length = [100, 100, 60, 100] #NEW LENGTH???????????!!!!!!!!!
 		#length = [100, 100, 100, 100]
		b = x[i] < length[i]	


		#Object+PSF Profile
		ax1 = fig.add_subplot(gs[0]) #since we only have 1 col, we have to only give the row into which the plot must go
		ax1.errorbar(x[i][b], object_prof[i][b], yerr=object_prof_unc[i][b], fmt='-o', color='darkblue') # markersize=0.5,
		ax1.errorbar(x[i][b], psf_prof[i][b], fmt='--', markersize=0.5, color='grey') #plot both obj and psf profiles in the set of ax1 axes PLOTTING	
		nbins = len(ax1.get_xticklabels())
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='lower'))
		#ax1.get_yaxis().set_label_coords(-0.01,0.5)
		plt.setp(ax1.get_xticklabels(), visible=False)
		ax1.set_yscale('log') #setting y axis to log scale
	
	
		#Residual Beam--> (Object - Beam)
		ax2 = fig.add_subplot(gs[1], sharex=ax1) 
		b[0:2]=False #We make the 0th and 1st point a false index since we don't want to plot the 1st two points in the res and fracres profiles. Always give the value one above the needed value for a range in [:].
		#b[0]=False #We make the 0th point a false index
		ax2.errorbar(x[i][b], res[i][b], yerr=object_prof_unc[i][b], fmt='--^', color='darkblue')			
		ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
		plt.setp(ax2.get_xticklabels(), visible=False)
		#nbins = len(ax1.get_xticklabels())
		ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max val on y axis to prevent merging with upper axis labels
		ax2.set_ylim(ymin=0) #starting y axis from zero
		#ax2.get_yaxis().set_label_coords(-0.1,0.5) 
		
	
		#Fractional Residual --> (Object - Beam)/Beam --> Residual as a fractional excess compared to Beam
		ax3 = fig.add_subplot(gs[2], sharex=ax1) 
		ax3.errorbar(x[i][b], fracres[i][b], yerr=fracres_unc[i][b], fmt='--s', color='darkblue')
		ax3.set_ylabel('Fractional Residual') #for log scaled y axis
		#ax3.get_yaxis().set_label_coords(-0.1,0.5)            
		ax3.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper')) 
		ax3.set_ylim(ymin=0)
		ax3.set_xlabel("Radius ($^{\prime\prime}$)")
		ax3.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
		
		
		
		#Labelling y axis of the Radial Profile and Residual Profile
		if i < 2: #the 0th and 1st 
			ax1.set_ylabel('Surface Brightness \n(Jy/arcsec$^2$)') #PACS in Jy/arcsec^2
			ax2.set_ylabel('Residual \n(Jy/arcsec$^2$)') #PACS in Jy/arcsec^2
		else:
			ax1.set_ylabel('Surface Brightness \n(mJy/arcsec$^2$)') #SCUBA2 in mJy/arcsec^2
			ax2.set_ylabel('Residual \n(mJy/arcsec$^2$)') #SCUBA2 in mJy/arcsec^2

		#if i == 1 or i == 3:
		#	ax1.yaxis.tick_right()
		#	ax1.yaxis.set_label_position("right")
		#	ax2.yaxis.tick_right()
		#	ax2.yaxis.set_label_position("right")
		#	ax3.yaxis.tick_right()
		#	ax3.yaxis.set_label_position("right")
		#invert (mirror image) y axis lable???????????? 
		#else: 
		#	ax1.yaxis.tick_left()
		#	ax1.yaxis.set_label_position("left")
		#	ax2.yaxis.tick_left()
		#	ax2.yaxis.set_label_position("left")
		#	ax3.yaxis.tick_left()
		#	ax3.yaxis.set_label_position("left")
			
		
	
		#Saving the figures #PUT THESE LINES BACK IN WHEN THE IMAGES NEED TO BE SAVED!!!!!!!!!!
		#save_file = star + str(wavelengths[i]) + '.pdf' #str() - converts i elememt of wavelengths in to a string. star is already a string. '.pdf' is a string. they all concatenate together with the +s.
		save_file = star + '_' + str(wavelengths[i]) + '.png'
		plt.savefig(save_file)
		
	#plt.show()
	
	


################ Saving Interpolated Residual data for SED plotting ############################################################

	
	
	#f = open(star+'_x_interp.dat', "w") #creating and wrtting(w) to a file. filename - given in mainbody. 
	#f.close()
	
	#print x_interp >> star+'_x_interp.dat'
	#f = open(star+'_x_interp.dat', "a") #append(a) the created data file 
	#for i in range(x_interp):
	#	outstring = str(x_interp[i]) + "\n"
	#	f.write(outstring)
	#f.close()
	

	
	#Resdidual Profile data for SED fitting for PACS data
	res_interp_SED = stellar_profile  - psf_profile 
	res_interp_unc_SED = np.sqrt((unc_profile**2)+((fractional_unc_psf*psf_profile)**2))
	

	#Resdidual Profile data for SED fitting where we change SCUBA2 values from mJy --> Jy	
	res_interp_SED[:,2:] = res_interp_SED[:,2:]/1000 # Here we change only the SCUBA2 array dimensions. :=change all data in y axis. 2:=but only change 2 and 3 (all from two onwards) in x axis. Here it's just two dimensions not really x and y directions. 
	res_interp_unc_SED[:,2:] = np.sqrt(((unc_profile[:,2:]/1000)**2)+((fractional_unc_psf[2:]*(psf_profile[:,2:]/1000))**2)) # As frac_unc_psf is only 1 line with 4 values, it's 1 D with four entries. So we just change all values from 2 (i.e: 2 and 3) for SCUBA2. 
		#fractional_unc_psf[i_lam] - devide by a thousand????????????????????CHECK!!!!!!!!!!!!!!!!!!!!!!

	#print res_interp_SED
	#print res_interp_unc_SED
	
	#Creating a .dat file and saving the data needed for SED fitting in it. 
	res_interp_save = res_interp_SED.reshape((len(res_interp_SED),4))
	res_interp_unc_save = res_interp_unc_SED.reshape((len(res_interp_unc_SED),4))
	#print res_interp_save
	#print res_interp_unc_save


	f = open(star+'_res_interp.dat', "w") #creating and wrtting(w) to a file. filename - given in mainbody. 
	f.close()
	f = open(star+'_res_interp.dat', "a") #append(a) the created data file 
	for i in range(res_interp_SED.shape[0]):
		outstring = "\t" + str(res_interp_save[i,0]) + "\t" + str(res_interp_save[i,1]) + "\t" + str(res_interp_save[i,2]) + "\t" + str(res_interp_save[i,3]) + "\t" + "\n" #"\t" = tab space. "\n" = start new line. str() = converting/casting a value into a string-Textual data in Python is handled with str objects. We're just printing the data as text format(allthough they are numbers) into the file.
		f.write(outstring)
	f.close()


	f = open(star+'_res_interp_unc.dat', "w") #creating and wrtting(w) to a file. filename - given in mainbody. 
	f.close()
	f = open(star+'_res_interp_unc.dat', "a") #append(a) the created data file 
	for i in range(res_interp_unc_SED.shape[0]):
		outstring = "\t" + str(res_interp_unc_save[i,0]) + "\t" + str(res_interp_unc_save[i,1]) + "\t" + str(res_interp_unc_save[i,2]) + "\t" + str(res_interp_unc_save[i,3]) + "\t" + "\n" #"\t" = tab space. "\n" = start new line. str() = converting/casting a value into a string-Textual data in Python is handled with str objects. We're just printing the data as text format(allthough they are numbers) into the file.
		f.write(outstring)
	f.close() 
	
   
   
        
########################### Plotting only profile and residual profile. #################################

	for i in range(len(x)):
		fig = plt.figure(figsize=(9,6)) #(width,height)
		gs = gridspec.GridSpec(2, 1, height_ratios=[0.5, 0.5]) #Note capital G and S in .GridSpec
		gs.update(hspace=0)
		length = [100, 100, 60, 100] #x axis arcseconds lengths for 70, 160, 450, 850 in this specific order!
		b = x[i] < length[i]	


		#Object+PSF Profile
		ax1 = fig.add_subplot(gs[0]) #since we only have 1 col, we have to only give the row into which the plot must go
		ax1.errorbar(x[i][b], object_prof[i][b], yerr=object_prof_unc[i][b], fmt='-o', color='darkblue') # markersize=0.5,
		ax1.errorbar(x[i][b], psf_prof[i][b], fmt='--', markersize=0.5, color='grey') #plot both obj and psf profiles in the set of ax1 axes PLOTTING	
		nbins = len(ax1.get_xticklabels())
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='lower'))
		#ax1.get_yaxis().set_label_coords(-0.01,0.5)
		plt.setp(ax1.get_xticklabels(), visible=False)
		ax1.set_yscale('log') #setting y axis to log scale
	
	
		#Residual Beam--> (Object - Beam)
		ax2 = fig.add_subplot(gs[1], sharex=ax1) 
		b[0:2]=False #We make the 0th and 1st point a false index since we don't want to plot the 1st two points in the res and fracres profiles. Always give the value one above the needed value for a range in [:].
		#b[0]=False #We make the 0th point a false index
		ax2.errorbar(x[i][b], res[i][b], yerr=object_prof_unc[i][b], fmt='--^', color='darkblue')			
		ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
		#plt.setp(ax2.get_xticklabels())
		#plt.setp(ax2.get_xticklabels(), visible=False)
		#nbins = len(ax1.get_xticklabels())
		ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max val on y axis to prevent merging with upper axis labels
		ax2.set_ylim(ymin=0) #starting y axis from zero
		ax2.set_xlabel("Radius ($^{\prime\prime}$)")
		
	
		#Fractional Residual --> (Object - Beam)/Beam --> Residual as a fractional excess compared to Beam
		#ax3 = fig.add_subplot(gs[2], sharex=ax1) 
		#ax3.errorbar(x[i][b], fracres[i][b], yerr=fracres_unc[i][b], fmt='--s', color='darkblue')
		#ax3.set_ylabel('Fractional Residual') #for log scaled y axis
		#ax3.get_yaxis().set_label_coords(-0.1,0.5)            
		#ax3.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper')) 
		#ax3.set_ylim(ymin=0)
		#ax3.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
		
		
		
		#Labelling y axis of the Radial Profile and Residual Profile
		if i < 2: #the 0th and 1st 
			ax1.set_ylabel('Surface Brightness \n(Jy/arcsec$^2$)') #PACS in Jy/arcsec^2
			ax2.set_ylabel('Residual \n(Jy/arcsec$^2$)') #PACS in Jy/arcsec^2
		else:
			ax1.set_ylabel('Surface Brightness \n(mJy/arcsec$^2$)') #SCUBA2 in mJy/arcsec^2
			ax2.set_ylabel('Residual \n(mJy/arcsec$^2$)') #SCUBA2 in mJy/arcsec^2

			
		
	
		#Saving the figures #PUT THESE LINES BACK IN WHEN THE IMAGES NEED TO BE SAVED!!!!!!!!!!
		#save_file = star + str(wavelengths[i]) + '.pdf' #str() - converts i elememt of wavelengths in to a string. star is already a string. '.pdf' is a string. they all concatenate together with the +s.
		save_file = star + '_' + str(wavelengths[i]) + '_onlyprofile+res.png'
		plt.savefig(save_file)
		plt.show()
	
	                
                
                
                
                
                
                
                
     






 
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
