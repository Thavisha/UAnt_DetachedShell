import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec
from scipy.integrate import simps, cumtrapz
from numpy import trapz
import matplotlib.ticker as mtick
from astropy import constants as const
from astropy import units as u
from astropy.analytic_functions import blackbody_lambda, blackbody_nu
from astropy.table import Table
from astropy.constants import c #c=speed of light 


"""
###################################################################################################################################
    ########### Dust and Gass Mass Loss ################################
###################################################################################################################################

- Calculate the total dust mass of the enhanced region of the density profile of U Ant (12" - 52") 

- Total dust mass derived by integrating over the density profile resulting from SED fitting. 
- projected size converted to physical size using terminal velocity and distance to source using the values presented in the paper. 
- Total gas mass derived using CO gas mass loss rates from DeBeck et al., 2010 and the age of CSE at the radial point. 
- All values derived up to the 3 sigma extension point of 850micron.  

- Monter Carlo Uncertainties Added to derive the Uncertainty on the total mass

###### Input files required ####################

1) 160_3sigma_Extensions_And_Total_Fluxes_And_DeBeckMassLossRates.csv - table containing 3 sigma extensions, DeBeck data and distances, etc needed. 
2) x_interp.csv - file containing x radial points
3) .csv files containing the density profiles derived from SED fitting section for all source listed in file 1). eg: cit6_density_radial_output.csv


##### Output ############

.csv table containing i)source, ii)3 sigma radius in arcsec, iii)3 sigma radius in cm, iv)3 sigma radius in pc, v)DeBeck mass loss rate, vi) v)DeBeck terminal velocity, vi)age of CSE at 3sigma radius, vii)total gas mass, viii)Dust mass using derived gas mass and 1/200 d:g ratio, ix) Total dust mass from integrating the density profile   

"""



#Total Mass from the Density Profile Derived from SEDs. #Integration of Density Profile
def total_mass_from_density_profile(star, three_sigma_cm, distance, nskip=3): #nskip=1 - Skip values in the 0th and 1st and 2nd radial point

	density_data_filename = star + '_density_radial_output.csv'
	density_data = Table.read(density_data_filename, format='ascii.csv')

	density_profile = density_data['Density']
	density_profile = density_profile[nskip:]
	density_profile_plusUnc = density_data['Plus_Unc']
	density_profile_plusUnc = density_profile_plusUnc[nskip:]
	density_profile_minusUnc = density_data['Minus_Unc']
	density_profile_minusUnc = density_profile_minusUnc[nskip:]

	density_profile_unc = (density_profile_plusUnc + density_profile_minusUnc)/2 #Added here to work for Monte Carlo. Why can't we use the plus and minus seperately when giving it to line 47 - - Because np.random.normal resamples on a symetric grid so you need symmetric uncertainties as given here. 

	x = density_data['x_interp'] #In arcsec
	x = x[nskip:]
	x_cm = (distance * x) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
	three_sigma_cm_limit = x_cm < three_sigma_cm



	#Adding Monter Carlo to get the Uncetainty on the Integrated Density
	Integrated_Density_list = [] #Creating an empty list to hold the Monte Carlo outputs

	for i in range(2000):

		#numpy.random.normal(loc=0.0, scale=1.0, size=None) # Loc =mean, scale = std.dev(=sigma), size = no of random sampling. In our case mean=flux, sigma=flux_unc. 
		#density_profile_resampled = np.random.normal(density_profile, [density_profile_plusUnc, density_profile_minusUnc]) Why can't we give the uncertianties this way - Because np.random.normal resamples on a symetric grid so you need symmetric uncertainties as given below. 
		density_profile_resampled = np.random.normal(density_profile, density_profile_unc)
		
		Integrated_Density = simps((2*np.pi*x_cm[three_sigma_cm_limit]*(10**density_profile_resampled[three_sigma_cm_limit])), x_cm[three_sigma_cm_limit]) #Units = g
		Integrated_Density = Integrated_Density /  1.9884754e33 #Convert to Solar Masses
		#print('Integrated_Density =', Integrated_Density) 
		#print (cumtrapz((2*np.pi*x_cm[three_sigma_cm_limit]*(10**density_profile[three_sigma_cm_limit])), x_cm[three_sigma_cm_limit])  /  1.9884754e33 )
		Integrated_Density_list.append(Integrated_Density)

	
	Final_Density_and_Unc = np.percentile(Integrated_Density_list, [16, 50, 84]) #Get the 16th(-unc), 50th(median=value), 84th(+unc) of the Periods

	Final_Density = Final_Density_and_Unc[1]

	Final_Density_Plus_Unc = Final_Density_and_Unc[2] - Final_Density_and_Unc[1] #+Unc = 84th Percentile - Median

	Final_Density_Minus_Unc = Final_Density_and_Unc[1] - Final_Density_and_Unc[0]  #-Unc = Median - 16th Percentile 

	Final_Unc = (Final_Density_Plus_Unc + Final_Density_Minus_Unc)/2

	print('Final_Period =', Final_Density, '+/-', Final_Unc, '(+', Final_Density_Plus_Unc, '-', Final_Density_Minus_Unc, ')') #Period

	
	return Final_Density, Final_Unc, Final_Density_Plus_Unc, Final_Density_Minus_Unc






if __name__=="__main__":

	#Heavily Edited for U Ant!! Don't use for anything else!!!!!!!!!!

	wavelength = '850' #Get masses for the 850 3 sigma cause you need at least two radial points to get an accurate density. so 70 extension is not enough. at 850 there's both 70 and 850 detections
	wavelength_val = 850
	data_filename = 'UAnt_Details.csv'
	input_data = Table.read(data_filename, format='ascii.csv')
	print(input_data.columns)
	calculated_densitites_output_filename = 'UAnt_TotalMass_output_withMonteCarloUncertainties.csv'
	#print source_list


	f = open(calculated_densitites_output_filename, 'w')
	f.write("source" + "," + wavelength+"_three_sigma(arcsec)" + "," + wavelength+"_three_sigma(cm)" + "," + wavelength+"_three_sigma(pc)" + "," + "DeBeck_etal_2010_Terminal_Velocity_fromCOdata(cm/s)" + "," + "3sigma_radius_850_Time(years)(Time=Dist/Vel)unit_sToYr=*3.171e-8" + "," + "Final_Intergrated_Dust_Mass_fromDensityProfile(Msun)" + "," + "Final_Uncetainty_on_Dust_Mass_fromMonteCarlo(Msun)" + "," + "Plus_Uncetainty_on_Dust_Mass_fromMonteCarlo(Msun)" + "," + "Minus_Uncetainty_on_Dust_Mass_fromMonteCarlo(Msun)" + "\n")
	f.close()

	for source in input_data: 

		star = str(source['Source'])
		distance = source['Distance(pc)'] #In pc
		distance_cm = distance * 3.086E+18 #Convert from pc to cm
		three_sigma_arcsec = source['three_sigma_Extension(arcsec)'] #In arcsec. needs to be converted to cm using the source distances
		three_sigma_cm = (distance * three_sigma_arcsec) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
		three_sigma_pc = three_sigma_cm * 3.24078E-019 
		#DeBeck_MassLoss_Rate =  source['DeBeck_CO_Mass_Loss_Rate(M_sun/yr)']
		DeBeck_Terminal_Velocity = source['DeBeck_etal_2010_Terminal_Velocity_fromCOdata(cm/s)']

		
		Final_Density, Final_Unc, Final_Density_Plus_Unc, Final_Density_Minus_Unc = total_mass_from_density_profile(star, three_sigma_cm, distance, nskip=1)

		three_sigma_radius_inTimeUnits = (three_sigma_cm/DeBeck_Terminal_Velocity)*3.171E-8 #3sigma_radius_850_Time (years) (Time=Dist/Vel) unit_s to Yr = *3.171e-8. Same as doing (seventy_3sigma_cm / velocity) / 3.154e+7 #/3.154e+7 to convert from s to years

		#DeBeck_total_CO_MassLoss_in_3sigma_time = DeBeck_MassLoss_Rate * three_sigma_radius_inTimeUnits

		#Dust_Mass_fromDeBeck = DeBeck_total_CO_MassLoss_in_3sigma_time * 0.005 #Accepted dust:gas ratio = 1/200 = 0.005

		#Dust_Mass_Ratio = Integrated_Density / Dust_Mass_fromDeBeck 

		#Dust_to_Gas_Ratio = Integrated_Density / DeBeck_total_CO_MassLoss_in_3sigma_time


		
		f = open(calculated_densitites_output_filename, 'a')
		outstring = str(star)  + ","  + str(three_sigma_arcsec)  + ","  + str(three_sigma_cm)   + ","  + str(three_sigma_pc) + "," + str(DeBeck_Terminal_Velocity) + "," + str(three_sigma_radius_inTimeUnits) + "," + str(Final_Density) + "," + str(Final_Unc) + "," + str(Final_Density_Plus_Unc) + "," + str(Final_Density_Minus_Unc) + "\n"
		f.write(outstring)	
		f.close()


		























