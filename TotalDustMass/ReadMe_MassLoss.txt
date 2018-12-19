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




