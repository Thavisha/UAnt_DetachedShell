import numpy as np
from astropy.table import Table
import warnings

"""
Calculating the Reduced Chisquared values of each model rad.prof at each wavelenght with respect to the observed rad prof.

When running this give an output filename in the command line to save the Chi Squared results in the file
eg: python RedChiSquared_Image.py > RedChiSquaredValues_Image.csv #This prints the output to the filename RedChiSquaredValues_Image.csv

"""


star = 'uant'
distance = 268.097  #Units = pc 

############# Observed Profiles ##################################################

#Loading and opening PACS 70 micron data files
Stellar_Data_70 = Table.read(star+'_UnInterpolated_70.csv', format='ascii.csv')
x_70 = Stellar_Data_70['x_axis(arcsec)']
#x_70_cm = (distance * x_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_70_pc = (distance * x_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_70 = Stellar_Data_70['Stellar_Profile(Jy/arcsec^2)']
stellar_70 = stellar_70 * 1000 #Converting from Jy to mJy
stellar_unc_70 = Stellar_Data_70['Stellar_Unc(Jy/arcsec^2)']
stellar_unc_70 = stellar_unc_70 * 1000 #Converting from Jy to mJy


#Loading and opening PACS 160 micron data files
Stellar_Data_160 = Table.read(star+'_UnInterpolated_160.csv', format='ascii.csv')
x_160 = Stellar_Data_160['x_axis(arcsec)']
x_160_cm = (distance * x_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
stellar_160 = Stellar_Data_160['Stellar_Profile(Jy/arcsec^2)']
stellar_160 = stellar_160 * 1000 #Converting from Jy to mJy
stellar_unc_160 = Stellar_Data_160['Stellar_Unc(Jy/arcsec^2)']
stellar_unc_160 = stellar_unc_160 * 1000 #Converting from Jy to mJy


#Loading and opening SCUBA2 450 micron data files
Stellar_Data_450 = Table.read(star+'_UnInterpolated_450.csv', format='ascii.csv')
x_450 = Stellar_Data_450['x_axis(arcsec)']
x_450_cm = (distance * x_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
stellar_450 = Stellar_Data_450['Stellar_Profile(mJy/arcsec^2)']
stellar_unc_450 = Stellar_Data_450['Stellar_Unc(mJy/arcsec^2)']


#Loading and opening SCUBA2 850 micron data files
Stellar_Data_850 = Table.read(star+'_UnInterpolated_850.csv', format='ascii.csv')
x_850 = Stellar_Data_850['x_axis(arcsec)']
x_850_cm = (distance * x_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
stellar_850 = Stellar_Data_850['Stellar_Profile(mJy/arcsec^2)']
stellar_unc_850 = Stellar_Data_850['Stellar_Unc(mJy/arcsec^2)']



############# Model Profiles for 70 micron (PACS BLUE) ##################################################

Stellar_Data_ContEm_Only_70 = Table.read('Model_ContEm_Only_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_ContEm_Only_70 = Stellar_Data_ContEm_Only_70['x_axis(arcsec)']
#x_ContEm_Only_70_cm = (distance * x_ContEm_Only_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_ContEm_Only_70_pc = (distance * x_ContEm_Only_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_ContEm_Only_70 = Stellar_Data_ContEm_Only_70['Stellar_Profile(Jy/arcsec^2)']
stellar_ContEm_Only_70 = stellar_ContEm_Only_70 * 1000 #Converting from Jy to mJy

Stellar_Data_InnerShell_70 = Table.read('Model_InnerShell_ContEm_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_InnerShell_70 = Stellar_Data_InnerShell_70['x_axis(arcsec)']
#x_InnerShell_70_cm = (distance * x_InnerShell_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_InnerShell_70_pc = (distance * x_InnerShell_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_InnerShell_70 = Stellar_Data_InnerShell_70['Stellar_Profile(Jy/arcsec^2)']
stellar_InnerShell_70 = stellar_InnerShell_70 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_K2010_70 = Table.read('Model_OuterShell_ContEm_K2010_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_OuterShell_K2010_70 = Stellar_Data_OuterShell_K2010_70['x_axis(arcsec)']
#x_OuterShell_K2010_70_cm = (distance * x_OuterShell_K2010_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_K2010_70_pc = (distance * x_OuterShell_K2010_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_K2010_70 = Stellar_Data_OuterShell_K2010_70['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_K2010_70 = stellar_OuterShell_K2010_70 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_M2010_70 = Table.read('Model_OuterShell_ContEm_M2010_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_OuterShell_M2010_70 = Stellar_Data_OuterShell_M2010_70['x_axis(arcsec)']
#x_OuterShell_M2010_70_cm = (distance * x_OuterShell_M2010_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_M2010_70_pc = (distance * x_OuterShell_M2010_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_M2010_70 = Stellar_Data_OuterShell_M2010_70['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_M2010_70 = stellar_OuterShell_M2010_70 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Even_70 = Table.read('Model_FourShells+ContEm_EvenMass_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_FourShells_Even_70 = Stellar_Data_FourShells_Even_70['x_axis(arcsec)']
#x_FourShells_Even_70_cm = (distance * x_FourShells_Even_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Even_70_pc = (distance * x_FourShells_Even_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Even_70 = Stellar_Data_FourShells_Even_70['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Even_70 = stellar_FourShells_Even_70 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Uneven_70 = Table.read('Model_FourShells_ContEm_UnevenMass_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_FourShells_Uneven_70 = Stellar_Data_FourShells_Uneven_70['x_axis(arcsec)']
#x_FourShells_Uneven_70_cm = (distance * x_FourShells_Uneven_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Uneven_70_pc = (distance * x_FourShells_Uneven_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Uneven_70 = Stellar_Data_FourShells_Uneven_70['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Uneven_70 = stellar_FourShells_Uneven_70 * 1000 #Converting from Jy to mJy

############# Model Profiles for 160 micron (PACS RED) ##################################################

Stellar_Data_ContEm_Only_160 = Table.read('Model_ContEm_Only_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_ContEm_Only_160 = Stellar_Data_ContEm_Only_160['x_axis(arcsec)']
#x_ContEm_Only_160_cm = (distance * x_ContEm_Only_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_ContEm_Only_160_pc = (distance * x_ContEm_Only_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_ContEm_Only_160 = Stellar_Data_ContEm_Only_160['Stellar_Profile(Jy/arcsec^2)']
stellar_ContEm_Only_160 = stellar_ContEm_Only_160 * 1000 #Converting from Jy to mJy

Stellar_Data_InnerShell_160 = Table.read('Model_InnerShell_ContEm_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_InnerShell_160 = Stellar_Data_InnerShell_160['x_axis(arcsec)']
#x_InnerShell_160_cm = (distance * x_InnerShell_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_InnerShell_160_pc = (distance * x_InnerShell_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_InnerShell_160 = Stellar_Data_InnerShell_160['Stellar_Profile(Jy/arcsec^2)']
stellar_InnerShell_160 = stellar_InnerShell_160 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_K2010_160 = Table.read('Model_OuterShell_ContEm_K2010_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_OuterShell_K2010_160 = Stellar_Data_OuterShell_K2010_160['x_axis(arcsec)']
#x_OuterShell_K2010_160_cm = (distance * x_OuterShell_K2010_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_K2010_160_pc = (distance * x_OuterShell_K2010_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_K2010_160 = Stellar_Data_OuterShell_K2010_160['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_K2010_160 = stellar_OuterShell_K2010_160 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_M2010_160 = Table.read('Model_OuterShell_ContEm_M2010_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_OuterShell_M2010_160 = Stellar_Data_OuterShell_M2010_160['x_axis(arcsec)']
#x_OuterShell_M2010_160_cm = (distance * x_OuterShell_M2010_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_M2010_160_pc = (distance * x_OuterShell_M2010_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_M2010_160 = Stellar_Data_OuterShell_M2010_160['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_M2010_160 = stellar_OuterShell_M2010_160 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Even_160 = Table.read('Model_FourShells+ContEm_EvenMass_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_FourShells_Even_160 = Stellar_Data_FourShells_Even_160['x_axis(arcsec)']
#x_FourShells_Even_160_cm = (distance * x_FourShells_Even_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Even_160_pc = (distance * x_FourShells_Even_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Even_160 = Stellar_Data_FourShells_Even_160['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Even_160 = stellar_FourShells_Even_160 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Uneven_160 = Table.read('Model_FourShells_ContEm_UnevenMass_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_FourShells_Uneven_160 = Stellar_Data_FourShells_Uneven_160['x_axis(arcsec)']
#x_FourShells_Uneven_160_cm = (distance * x_FourShells_Uneven_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Uneven_160_pc = (distance * x_FourShells_Uneven_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Uneven_160 = Stellar_Data_FourShells_Uneven_160['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Uneven_160 = stellar_FourShells_Uneven_160 * 1000 #Converting from Jy to mJy

############# Model Profiles for 450 micron (JCMT SCUBA 2 450) ##################################################

Stellar_Data_ContEm_Only_450 = Table.read('Model_ContEm_Only_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_ContEm_Only_450 = Stellar_Data_ContEm_Only_450['x_axis(arcsec)']
#x_ContEm_Only_450_cm = (distance * x_ContEm_Only_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_ContEm_Only_450_pc = (distance * x_ContEm_Only_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_ContEm_Only_450 = Stellar_Data_ContEm_Only_450['Stellar_Profile(Jy/arcsec^2)']
stellar_ContEm_Only_450 = stellar_ContEm_Only_450 * 1000 #Converting from Jy to mJy

Stellar_Data_InnerShell_450 = Table.read('Model_InnerShell_ContEm_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_InnerShell_450 = Stellar_Data_InnerShell_450['x_axis(arcsec)']
#x_InnerShell_450_cm = (distance * x_InnerShell_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_InnerShell_450_pc = (distance * x_InnerShell_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_InnerShell_450 = Stellar_Data_InnerShell_450['Stellar_Profile(Jy/arcsec^2)']
stellar_InnerShell_450 = stellar_InnerShell_450 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_K2010_450 = Table.read('Model_OuterShell_ContEm_K2010_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_OuterShell_K2010_450 = Stellar_Data_OuterShell_K2010_450['x_axis(arcsec)']
#x_OuterShell_K2010_450_cm = (distance * x_OuterShell_K2010_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_K2010_450_pc = (distance * x_OuterShell_K2010_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_K2010_450 = Stellar_Data_OuterShell_K2010_450['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_K2010_450 = stellar_OuterShell_K2010_450 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_M2010_450 = Table.read('Model_OuterShell_ContEm_M2010_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_OuterShell_M2010_450 = Stellar_Data_OuterShell_M2010_450['x_axis(arcsec)']
#x_OuterShell_M2010_450_cm = (distance * x_OuterShell_M2010_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_M2010_450_pc = (distance * x_OuterShell_M2010_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_M2010_450 = Stellar_Data_OuterShell_M2010_450['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_M2010_450 = stellar_OuterShell_M2010_450 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Even_450 = Table.read('Model_FourShells+ContEm_EvenMass_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_FourShells_Even_450 = Stellar_Data_FourShells_Even_450['x_axis(arcsec)']
#x_FourShells_Even_450_cm = (distance * x_FourShells_Even_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Even_450_pc = (distance * x_FourShells_Even_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Even_450 = Stellar_Data_FourShells_Even_450['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Even_450 = stellar_FourShells_Even_450 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Uneven_450 = Table.read('Model_FourShells_ContEm_UnevenMass_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_FourShells_Uneven_450 = Stellar_Data_FourShells_Uneven_450['x_axis(arcsec)']
#x_FourShells_Uneven_450_cm = (distance * x_FourShells_Uneven_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Uneven_450_pc = (distance * x_FourShells_Uneven_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Uneven_450 = Stellar_Data_FourShells_Uneven_450['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Uneven_450 = stellar_FourShells_Uneven_450 * 1000 #Converting from Jy to mJy

############# Model Profiles for 850 micron (JCMT SCUBA 2 850) ##################################################

Stellar_Data_ContEm_Only_850 = Table.read('Model_ContEm_Only_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_ContEm_Only_850 = Stellar_Data_ContEm_Only_850['x_axis(arcsec)']
#x_ContEm_Only_850_cm = (distance * x_ContEm_Only_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_ContEm_Only_850_pc = (distance * x_ContEm_Only_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_ContEm_Only_850 = Stellar_Data_ContEm_Only_850['Stellar_Profile(Jy/arcsec^2)']
stellar_ContEm_Only_850 = stellar_ContEm_Only_850 * 1000 #Converting from Jy to mJy

Stellar_Data_InnerShell_850 = Table.read('Model_InnerShell_ContEm_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_InnerShell_850 = Stellar_Data_InnerShell_850['x_axis(arcsec)']
#x_InnerShell_850_cm = (distance * x_InnerShell_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_InnerShell_850_pc = (distance * x_InnerShell_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_InnerShell_850 = Stellar_Data_InnerShell_850['Stellar_Profile(Jy/arcsec^2)']
stellar_InnerShell_850 = stellar_InnerShell_850 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_K2010_850 = Table.read('Model_OuterShell_ContEm_K2010_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_OuterShell_K2010_850 = Stellar_Data_OuterShell_K2010_850['x_axis(arcsec)']
#x_OuterShell_K2010_850_cm = (distance * x_OuterShell_K2010_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_K2010_850_pc = (distance * x_OuterShell_K2010_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_K2010_850 = Stellar_Data_OuterShell_K2010_850['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_K2010_850 = stellar_OuterShell_K2010_850 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_M2010_850 = Table.read('Model_OuterShell_ContEm_M2010_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_OuterShell_M2010_850 = Stellar_Data_OuterShell_M2010_850['x_axis(arcsec)']
#x_OuterShell_M2010_850_cm = (distance * x_OuterShell_M2010_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_M2010_850_pc = (distance * x_OuterShell_M2010_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_M2010_850 = Stellar_Data_OuterShell_M2010_850['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_M2010_850 = stellar_OuterShell_M2010_850 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Even_850 = Table.read('Model_FourShells+ContEm_EvenMass_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_FourShells_Even_850 = Stellar_Data_FourShells_Even_850['x_axis(arcsec)']
#x_FourShells_Even_850_cm = (distance * x_FourShells_Even_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Even_850_pc = (distance * x_FourShells_Even_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Even_850 = Stellar_Data_FourShells_Even_850['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Even_850 = stellar_FourShells_Even_850 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Uneven_850 = Table.read('Model_FourShells_ContEm_UnevenMass_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_FourShells_Uneven_850 = Stellar_Data_FourShells_Uneven_850['x_axis(arcsec)']
#x_FourShells_Uneven_850_cm = (distance * x_FourShells_Uneven_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Uneven_850_pc = (distance * x_FourShells_Uneven_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Uneven_850 = Stellar_Data_FourShells_Uneven_850['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Uneven_850 = stellar_FourShells_Uneven_850 * 1000 #Converting from Jy to mJy



##### Reduced ChiSqaured Generation #################

print("Wavelength and Model Type, Reduced ChiSqaured Result") #Prints a heading in the file given when running python. 


#Reduced RedChiSqaured for 70 micron 
radial_limit_obs_70 = x_70 < 100
radial_limit_Mod_70 = x_ContEm_Only_70 < 100

#Model Outer Shell (M2010) + Cont Emission - Shell size from Maercker+2010
OuterShell_M2010_70_RedChiSq = ( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_OuterShell_M2010_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) ) 
                                                                                                                                                / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
#print(np.shape(stellar_70[radial_limit_obs_70]))
print("outer_M2010_RedChiSq_70,", OuterShell_M2010_70_RedChiSq) 

#Model Outer Shell (K2010) + Cont Emission - Shell size from Kerschbaum+2010
OuterShell_K2010_70_RedChiSq =( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_OuterShell_K2010_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) )
                                                                                                                                    / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
print("outer_K2010_RedChiSq_70,", OuterShell_K2010_70_RedChiSq) 

#Model Inner Shell + Cont Emission
InnerShell_70_RedChiSq = ( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_InnerShell_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) )
                                                                                                                                / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
print("Inner_RedChiSq_70,", InnerShell_70_RedChiSq) 

#Model Four Shells + Cont Emission - Dust Mass Evenly spread between four shells
FourShells_Even_70_RedChiSq =( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_FourShells_Even_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) )
                                                                                                                                        / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
print("FourShell_Even_RedChiSq_70,", FourShells_Even_70_RedChiSq) 

#Model Four Shells + Cont Emission - Dust Mass Unevenly spread between four shells
FourShells_Uneven_70_RedChiSq = ( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_FourShells_Uneven_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) )
                                                                                                                                    / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
print("FourShell_Uneven_RedChiSq_70,", FourShells_Uneven_70_RedChiSq) 

#Model No Shells - Only Cont emission
ContEm_Only_70_RedChiSq = ( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_ContEm_Only_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) )
                                                                                                                                                    / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
print("ContEm_Only_RedChiSq_70,", ContEm_Only_70_RedChiSq ) 



#Reduced RedChiSqaured for 160 micron 
radial_limit_obs_160 = x_160 < 100
radial_limit_Mod_160 = x_ContEm_Only_160 < 100

#Model Outer Shell (M2010) + Cont Emission - Shell size from Maercker+2010
OuterShell_M2010_160_RedChiSq = ( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_OuterShell_M2010_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) )
                                                                                                                            / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
print("outer_M2010_RedChiSq_160,", OuterShell_M2010_160_RedChiSq) 

#Model Outer Shell (K2010) + Cont Emission - Shell size from Kerschbaum+2010
OuterShell_K2010_160_RedChiSq = ( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_OuterShell_K2010_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) )
                                                                                                                                / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
print("outer_K2010_RedChiSq_160,", OuterShell_K2010_160_RedChiSq) 

#Model Inner Shell + Cont Emission
InnerShell_160_RedChiSq = ( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_InnerShell_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) )
                                                                                                                                         / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
print("Inner_RedChiSq_160,", InnerShell_160_RedChiSq) 

#Model Four Shells + Cont Emission - Dust Mass Evenly spread between four shells
FourShells_Even_160_RedChiSq = ( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_FourShells_Even_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) )
                                                                                                                                             / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
print("FourShell_Even_RedChiSq_160,", FourShells_Even_160_RedChiSq) 

#Model Four Shells + Cont Emission - Dust Mass Unevenly spread between four shells
FourShells_Uneven_160_RedChiSq = ( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_FourShells_Uneven_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) )
                                                                                                                                      / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
print("FourShell_Uneven_RedChiSq_160,", FourShells_Uneven_160_RedChiSq) 

#Model No Shells - Only Cont emission
ContEm_Only_160_RedChiSq = ( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_ContEm_Only_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) )
                                                                                                                             / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
print("Only_RedChiSq_160,", ContEm_Only_160_RedChiSq ) 



#Reduced RedChiSqaured for 450 micron 
radial_limit_obs_450 = x_450 < 60
radial_limit_Mod_450 = x_ContEm_Only_450 < 60

#Model Outer Shell (M2010) + Cont Emission - Shell size from Maercker+2010
OuterShell_M2010_450_RedChiSq = ( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_OuterShell_M2010_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) )
                                                                                                                                                / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) )
print("outer_M2010_RedChiSq_450,", OuterShell_M2010_450_RedChiSq) 

#Model Outer Shell (K2010) + Cont Emission - Shell size from Kerschbaum+2010
OuterShell_K2010_450_RedChiSq = ( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_OuterShell_K2010_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) )
                                                                                                                    / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) )                                     
print("outer_K2010_RedChiSq_450,", OuterShell_K2010_450_RedChiSq) 

#Model Inner Shell + Cont Emission
InnerShell_450_RedChiSq = ( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_InnerShell_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) )
                                                                                                                / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) )
print("Inner_RedChiSq_450,", InnerShell_450_RedChiSq) 

#Model Four Shells + Cont Emission - Dust Mass Evenly spread between four shells
FourShells_Even_450_RedChiSq = ( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_FourShells_Even_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) )
                                                                                                                            / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) )
print("FourShell_Even_RedChiSq_450,", FourShells_Even_450_RedChiSq) 

#Model Four Shells + Cont Emission - Dust Mass Unevenly spread between four shells
FourShells_Uneven_450_RedChiSq = ( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_FourShells_Uneven_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) )
                                                                                                                                                     / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) ) 
print("FourShell_Uneven_RedChiSq_450,", FourShells_Uneven_450_RedChiSq) 

#Model No Shells - Only Cont emission
ContEm_Only_450_RedChiSq = ( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_ContEm_Only_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) )
                                                                                                                / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) )
print("ContEm_Only_RedChiSq_450,", ContEm_Only_450_RedChiSq ) 



#Reduced RedChiSqaured for 850 micron 
radial_limit_obs_850 = x_850 < 100
radial_limit_Mod_850 = x_ContEm_Only_850 < 100

#Model Outer Shell (M2010) + Cont Emission - Shell size from Maercker+2010
OuterShell_M2010_850_RedChiSq = ( ( np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_OuterShell_M2010_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 ) )
                                                                                                                            / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
print("outer_M2010_RedChiSq_850,", OuterShell_M2010_850_RedChiSq) 

#Model Outer Shell (K2010) + Cont Emission - Shell size from Kerschbaum+2010
OuterShell_K2010_850_RedChiSq = ( ( np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_OuterShell_K2010_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 )) 
                                                                                                                        / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
print("outer_K2010_RedChiSq_850,", OuterShell_K2010_850_RedChiSq) 

#Model Inner Shell + Cont Emission
InnerShell_850_RedChiSq = ( (np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_InnerShell_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 ))
                                                                                                                        / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
print("Inner_RedChiSq_850,", InnerShell_850_RedChiSq) 

#Model Four Shells + Cont Emission - Dust Mass Evenly spread between four shells
FourShells_Even_850_RedChiSq = ( (np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_FourShells_Even_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 ))
                                                                                                            / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
print("FourShell_Even_RedChiSq_850,", FourShells_Even_850_RedChiSq) 

#Model Four Shells + Cont Emission - Dust Mass Unevenly spread between four shells
FourShells_Uneven_850_RedChiSq = ( (np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_FourShells_Uneven_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 ))
                                                                                                                                        / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
print("FourShell_Uneven_RedChiSq_850,", FourShells_Uneven_850_RedChiSq) 

#Model No Shells - Only Cont emission
ContEm_Only_850_RedChiSq = ( (np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_ContEm_Only_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 ))
                                                                                                                                    / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
print("ContEm_Only_RedChiSq_850,", ContEm_Only_850_RedChiSq ) 














































