import matplotlib.pyplot as plt
import numpy as np
from numpy import trapz
from scipy.integrate import simps
from astropy.table import Table

"""
###################################################################################################################################
    ########### Uniform Mass Loss model ################################
###################################################################################################################################

- Generating a radially dependent dust mass column density profile assuming uniform mass loss. 

- We assume dust mass column density proportional to 1/r^2 relationship and extending it in the z and impact parameter directions. 


###### Input files required ####################

 - x_interp.dat file containing x axis information for upto which to generate the UMLM

###### Output ##################################  

1) .csv file containing the dust mass column density radial profile derived assuming uniform mass loss. This profile can now be over-plotted with the dust mass column density derived using SED fitting to show the deviation of the SED derived profile from uniform mass loss. 

"""

def Uniform_MassLoss_Model(impact_para, impact_para_0, z_0, r_max): 
#Project radius vs. Density is what we get out. Using density is proportional to 1/r^2 law and extending it in the z and impact parameter directions,

	impact_para = impact_para

	z_max = np.sqrt((r_max ** 2) - (impact_para ** 2))

	z = np.arange(0,z_max)

	density = (
			(np.sqrt((impact_para_0 ** 2) + (z_0 ** 2)))

				/(np.sqrt((impact_para ** 2) + (z ** 2)))

			) ** 2

	return density, z

f = open('Uniform_MassLoss_Model.csv', "w") #All other sources except IRC+10216
#f = open('Uniform_MassLoss_Model_longXaxisforIRC10216.csv', "w") #Long x axis for IRC+10216
f.write("x_interp" + "," + "Uniform_MassLoss_Model_Values" + "\n") 
f.close()



if __name__=="__main__":

	surface_density = [] #Defining an empty list to save the surface density output

	x_interp_data = np.loadtxt('x_interp_Uniform_MassLoss_Model.dat') # All sources except IRC+10216
	#x_interp_data = np.loadtxt('x_interp_irc10216_longxaxis_Uniform_MassLoss_Model.dat') # Long x axis for IRC+10216
	x_interp = x_interp_data.reshape([len(x_interp_data),1]) 

	for value in x_interp[:-1]: 
		impact_para = value

		impact_para_0 = 20
		z_0 = impact_para_0

		r_max = x_interp[-1]

		density, z = Uniform_MassLoss_Model(impact_para, impact_para_0, z_0, r_max)				
		
		surface_density.append(simps(density, z)) #Calculating surface density and saving it onto the blank surface_density list
	
	f = open('Uniform_MassLoss_Model.csv', "a") #All other sources except IRC+10216
	#f = open('Uniform_MassLoss_Model_longXaxisforIRC10216.csv', "a") #Long x axis for IRC+10216
	for i in range(len(surface_density)):
		outstring = str(x_interp[i,0]) + "," + str(surface_density[i]) + "\n" 
		f.write(outstring)
	f.close()
	
	#print x_interp[:,0].shape
	#print x_interp[:,0], surface_density,x_interp[1,0]

	



	




