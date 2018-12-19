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
Generating image of the collpased HARP image of U Ant
"""




Source = Table.read('UAnt_ContourPlot_Info.csv', format='ascii.csv')
fig = aplpy.FITSFigure('UAnt_HARP_Collapsed.fits')

star = 'UAnt'
x_center = 10 #Units in pixels 
y_center = 10  #Units in pixels 
beam_fwhm_harp = 14/3600.


cmap_type = 'cubehelix'
vmin_val = -2.2812524
vmax_val = 10
#vmid_val = 0.00776864
stretch_type = 'sqrt'
fig_size = 56/3600. #||150" in degrees

#Position to put beam shape at the edge
x_corner_pix = 2 #Units in pixels 
y_corner_pix = 2  #Units in pixels 

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
fig.colorbar.set_axis_label_text('K km s$^{-1}$')
fig.colorbar.set_axis_label_pad(25)
fig.colorbar.set_axis_label_rotation(270)
#fig.colorbar.show()

#Plotting beam on image
beam_radius = beam_fwhm_harp/2 #rad = FWHM/2 (diam=fwhm) #units=degrees (radius 13"=0.00361degrees). #For Radius divide by 2!!!!!!!!!!!!!! FWHM=Daimeter!


#Plotting beam as a circle at the edge of the image
fig.show_circles(x_corner, y_corner, beam_radius, edgecolor='none', facecolor='black', dashes='--', linewidth=2) #Prim.Beam (xcenter, ycenter, radius) 




fig.save('UAnt_HARP_ContourPlot.png', dpi=300)
plt.show()







