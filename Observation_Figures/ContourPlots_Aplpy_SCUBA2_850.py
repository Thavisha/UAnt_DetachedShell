import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
import aplpy
from astropy import units as u
from astropy.table import Table

font = {'family' : 'normal',
        'size'   : 16,
	'weight' : 'medium'}

plt.rc('font', **font)

"""
Generating image of the SCUBA-2 850micron image of U Ant

"""



wavelength = 850 
Source = Table.read('UAnt_ContourPlot_Info_850.csv', format='ascii.csv')
fig = aplpy.FITSFigure('uant_850.fits')

star = Source['Source']
x_center = Source['Center_X'] #Units in pixels 
y_center = Source['Center_Y']  #Units in pixels 
three_sigmaRadius = Source['ThreeSigma_Radius(arcsec)'] #Units =arcsec
Shell_1 = Source['Shell1_arcsec'] #Units =arcsec
Shell_2 = Source['Shell2_arcsec'] #Units =arcsec
Shell_3 = Source['Shell3_arcsec'] #Units =arcsec


beam_fwhm_dict = {'850': 13/3600., '450': 7.9/3600.}
beam_fwhm = beam_fwhm_dict[str(wavelength)] #0.00219444 #SCUBA2 Beam FWHM. #Units=degrees #850micron main beam =13"=0.00361degrees #450micron main beam =7.9"=0.00219444degrees


cmap_type = 'cubehelix'
vmin_val = 0
vmax_val = 0.1633763
#vmid_val = 0.00776864
stretch_type = 'sqrt'
fig_size = 150/3600. #||150" in degrees

#Position to put beam shape at the edge
x_corner_pix = 100  #Units in pixels 
y_corner_pix = 90  #Units in pixels 

#Creating Plot using Aplpy

fig.show_colorscale(vmin=vmin_val, vmax=vmax_val, cmap=cmap_type, stretch=stretch_type) #vmid=vmid_val,



#Cropping plot to the size we want AND Converting central position to wcs coord	
x_world, y_world = fig.pixel2world(x_center, y_center) #Converting pixels into degrees for plotting 
fig.recenter(x_world, y_world, fig_size) #Units should be in degrees!! || (x, y, size_in_degrees) 


#Converting beam position to wcs coord
x_corner, y_corner = fig.pixel2world(x_corner_pix, y_corner_pix)


#Adding Colourbar
fig.add_colorbar()
#fig.colorbar.set_width(0.2)
#fig.colorbar.set_pad(0.5)
#fig.colorbar.set_location('right')
fig.colorbar.set_axis_label_text('Surface Brightness (mJy arcsec$^{-2}$)')
fig.colorbar.set_axis_label_pad(25)
fig.colorbar.set_axis_label_rotation(270)
#fig.colorbar.show()

#Plotting beam on image
beam_radius = beam_fwhm/2 #rad = FWHM/2 (diam=fwhm) #units=degrees (radius 13"=0.00361degrees). #For Radius divide by 2!!!!!!!!!!!!!! FWHM=Daimeter!
three_sigmaRadius_degrees = three_sigmaRadius / 3600. #Converting from arcseconds to degrees.  
shell1 = Shell_1 / 3600.
shell2 = Shell_2 / 3600.
shell3 = Shell_3 / 3600.

#Plotting beam as a circle at the edge of the image
fig.show_circles(x_corner, y_corner, beam_radius, edgecolor='none', facecolor='white', dashes='--', linewidth=2) #Prim.Beam (xcenter, ycenter, radius) 

#Plotting 3sigma Radius as a circle
fig.show_circles(x_world, y_world, three_sigmaRadius_degrees, edgecolor='white', dashes='--', linewidth=2) #3sig.Radius (xcenter, ycenter, radius) 

#Optical Shells
#fig.show_circles(x_world, y_world, shell1, edgecolor='white', dashes='--', linewidth=2) #Shell 1 24" (xcenter, ycenter, radius)
#fig.show_circles(x_world, y_world, shell2, edgecolor='white', dashes=':', linewidth=2) #Shell 2 37" (xcenter, ycenter, radius) 
#fig.show_circles(x_world, y_world, shell3, edgecolor='white', dashes=':', linewidth=2) #Shell 2 43" (xcenter, ycenter, radius) 


fig.save('UAnt_850_ContourPlot.png', dpi=300)
plt.show()







