import numpy as np
from astropy.table import Table
import warnings

"""
Calculating the Reduced Chisquared values of each model rad.prof at each wavelength with respect to the observed rad prof.

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

Stellar_Data_InnerShell_70 = Table.read('Model_InnerShell_Only_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_InnerShell_70 = Stellar_Data_InnerShell_70['x_axis(arcsec)']
#x_InnerShell_70_cm = (distance * x_InnerShell_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_InnerShell_70_pc = (distance * x_InnerShell_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_InnerShell_70 = Stellar_Data_InnerShell_70['Stellar_Profile(Jy/arcsec^2)']
stellar_InnerShell_70 = stellar_InnerShell_70 * 1000 #Converting from Jy to mJy

Stellar_Data_ShellTwo_70 = Table.read('Model_ShellTwo_Only_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_ShellTwo_70 = Stellar_Data_ShellTwo_70['x_axis(arcsec)']
#x_ShellTwo_70_cm = (distance * x_ShellTwo_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_ShellTwo_70_pc = (distance * x_ShellTwo_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_ShellTwo_70 = Stellar_Data_ShellTwo_70['Stellar_Profile(Jy/arcsec^2)']
stellar_ShellTwo_70 = stellar_ShellTwo_70 * 1000 #Converting from Jy to mJy

Stellar_Data_ShellThree_70 = Table.read('Model_ShellThree_Only_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_ShellThree_70 = Stellar_Data_ShellThree_70['x_axis(arcsec)']
#x_ShellThree_70_cm = (distance * x_ShellThree_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_ShellThree_70_pc = (distance * x_ShellThree_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_ShellThree_70 = Stellar_Data_ShellThree_70['Stellar_Profile(Jy/arcsec^2)']
stellar_ShellThree_70 = stellar_ShellThree_70 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_M2010_70 = Table.read('Model_OuterShell_Only_M2010_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_OuterShell_M2010_70 = Stellar_Data_OuterShell_M2010_70['x_axis(arcsec)']
#x_OuterShell_M2010_70_cm = (distance * x_OuterShell_M2010_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_M2010_70_pc = (distance * x_OuterShell_M2010_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_M2010_70 = Stellar_Data_OuterShell_M2010_70['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_M2010_70 = stellar_OuterShell_M2010_70 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_K2010_70 = Table.read('Model_OuterShell_Only_K2010_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_OuterShell_K2010_70 = Stellar_Data_OuterShell_K2010_70['x_axis(arcsec)']
#x_OuterShell_K2010_70_cm = (distance * x_OuterShell_K2010_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_K2010_70_pc = (distance * x_OuterShell_K2010_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_K2010_70 = Stellar_Data_OuterShell_K2010_70['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_K2010_70 = stellar_OuterShell_K2010_70 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Uneven_70 = Table.read('Model_FourShells_UnEvenMass_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_FourShells_Uneven_70 = Stellar_Data_FourShells_Uneven_70['x_axis(arcsec)']
#x_FourShells_Uneven_70_cm = (distance * x_FourShells_Uneven_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Uneven_70_pc = (distance * x_FourShells_Uneven_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Uneven_70 = Stellar_Data_FourShells_Uneven_70['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Uneven_70 = stellar_FourShells_Uneven_70 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Even_70 = Table.read('Model_FourShells_EvenMass_HERSCHEL_PACS_BLUE_RadialProfile.csv', format='ascii.csv')
x_FourShells_Even_70 = Stellar_Data_FourShells_Even_70['x_axis(arcsec)']
#x_FourShells_Even_70_cm = (distance * x_FourShells_Even_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Even_70_pc = (distance * x_FourShells_Even_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Even_70 = Stellar_Data_FourShells_Even_70['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Even_70 = stellar_FourShells_Even_70 * 1000 #Converting from Jy to mJy



############# Model Profiles for 160 micron (PACS RED) ##################################################

Stellar_Data_InnerShell_160 = Table.read('Model_InnerShell_Only_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_InnerShell_160 = Stellar_Data_InnerShell_160['x_axis(arcsec)']
#x_InnerShell_160_cm = (distance * x_InnerShell_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_InnerShell_160_pc = (distance * x_InnerShell_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_InnerShell_160 = Stellar_Data_InnerShell_160['Stellar_Profile(Jy/arcsec^2)']
stellar_InnerShell_160 = stellar_InnerShell_160 * 1000 #Converting from Jy to mJy

Stellar_Data_ShellTwo_160 = Table.read('Model_ShellTwo_Only_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_ShellTwo_160 = Stellar_Data_ShellTwo_160['x_axis(arcsec)']
#x_ShellTwo_160_cm = (distance * x_ShellTwo_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_ShellTwo_160_pc = (distance * x_ShellTwo_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_ShellTwo_160 = Stellar_Data_ShellTwo_160['Stellar_Profile(Jy/arcsec^2)']
stellar_ShellTwo_160 = stellar_ShellTwo_160 * 1000 #Converting from Jy to mJy

Stellar_Data_ShellThree_160 = Table.read('Model_ShellThree_Only_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_ShellThree_160 = Stellar_Data_ShellThree_160['x_axis(arcsec)']
#x_ShellThree_160_cm = (distance * x_ShellThree_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_ShellThree_160_pc = (distance * x_ShellThree_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_ShellThree_160 = Stellar_Data_ShellThree_160['Stellar_Profile(Jy/arcsec^2)']
stellar_ShellThree_160 = stellar_ShellThree_160 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_M2010_160 = Table.read('Model_OuterShell_Only_M2010_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_OuterShell_M2010_160 = Stellar_Data_OuterShell_M2010_160['x_axis(arcsec)']
#x_OuterShell_M2010_160_cm = (distance * x_OuterShell_M2010_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_M2010_160_pc = (distance * x_OuterShell_M2010_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_M2010_160 = Stellar_Data_OuterShell_M2010_160['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_M2010_160 = stellar_OuterShell_M2010_160 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_K2010_160 = Table.read('Model_OuterShell_Only_K2010_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_OuterShell_K2010_160 = Stellar_Data_OuterShell_K2010_160['x_axis(arcsec)']
#x_OuterShell_K2010_160_cm = (distance * x_OuterShell_K2010_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_K2010_160_pc = (distance * x_OuterShell_K2010_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_K2010_160 = Stellar_Data_OuterShell_K2010_160['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_K2010_160 = stellar_OuterShell_K2010_160 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Uneven_160 = Table.read('Model_FourShells_UnEvenMass_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_FourShells_Uneven_160 = Stellar_Data_FourShells_Uneven_160['x_axis(arcsec)']
#x_FourShells_Uneven_160_cm = (distance * x_FourShells_Uneven_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Uneven_160_pc = (distance * x_FourShells_Uneven_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Uneven_160 = Stellar_Data_FourShells_Uneven_160['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Uneven_160 = stellar_FourShells_Uneven_160 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Even_160 = Table.read('Model_FourShells_EvenMass_HERSCHEL_PACS_RED_RadialProfile.csv', format='ascii.csv')
x_FourShells_Even_160 = Stellar_Data_FourShells_Even_160['x_axis(arcsec)']
#x_FourShells_Even_160_cm = (distance * x_FourShells_Even_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Even_160_pc = (distance * x_FourShells_Even_160) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Even_160 = Stellar_Data_FourShells_Even_160['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Even_160 = stellar_FourShells_Even_160 * 1000 #Converting from Jy to mJy

############# Model Profiles for 450 micron (JCMT SCUBA 2 450) ##################################################

Stellar_Data_InnerShell_450 = Table.read('Model_InnerShell_Only_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_InnerShell_450 = Stellar_Data_InnerShell_450['x_axis(arcsec)']
#x_InnerShell_450_cm = (distance * x_InnerShell_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_InnerShell_450_pc = (distance * x_InnerShell_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_InnerShell_450 = Stellar_Data_InnerShell_450['Stellar_Profile(Jy/arcsec^2)']
stellar_InnerShell_450 = stellar_InnerShell_450 * 1000 #Converting from Jy to mJy

Stellar_Data_ShellTwo_450 = Table.read('Model_ShellTwo_Only_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_ShellTwo_450 = Stellar_Data_ShellTwo_450['x_axis(arcsec)']
#x_ShellTwo_450_cm = (distance * x_ShellTwo_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_ShellTwo_450_pc = (distance * x_ShellTwo_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_ShellTwo_450 = Stellar_Data_ShellTwo_450['Stellar_Profile(Jy/arcsec^2)']
stellar_ShellTwo_450 = stellar_ShellTwo_450 * 1000 #Converting from Jy to mJy

Stellar_Data_ShellThree_450 = Table.read('Model_ShellThree_Only_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_ShellThree_450 = Stellar_Data_ShellThree_450['x_axis(arcsec)']
#x_ShellThree_450_cm = (distance * x_ShellThree_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_ShellThree_450_pc = (distance * x_ShellThree_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_ShellThree_450 = Stellar_Data_ShellThree_450['Stellar_Profile(Jy/arcsec^2)']
stellar_ShellThree_450 = stellar_ShellThree_450 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_M2010_450 = Table.read('Model_OuterShell_Only_M2010_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_OuterShell_M2010_450 = Stellar_Data_OuterShell_M2010_450['x_axis(arcsec)']
#x_OuterShell_M2010_450_cm = (distance * x_OuterShell_M2010_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_M2010_450_pc = (distance * x_OuterShell_M2010_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_M2010_450 = Stellar_Data_OuterShell_M2010_450['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_M2010_450 = stellar_OuterShell_M2010_450 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_K2010_450 = Table.read('Model_OuterShell_Only_K2010_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_OuterShell_K2010_450 = Stellar_Data_OuterShell_K2010_450['x_axis(arcsec)']
#x_OuterShell_K2010_450_cm = (distance * x_OuterShell_K2010_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_K2010_450_pc = (distance * x_OuterShell_K2010_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_K2010_450 = Stellar_Data_OuterShell_K2010_450['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_K2010_450 = stellar_OuterShell_K2010_450 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Uneven_450 = Table.read('Model_FourShells_UnEvenMass_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_FourShells_Uneven_450 = Stellar_Data_FourShells_Uneven_450['x_axis(arcsec)']
#x_FourShells_Uneven_450_cm = (distance * x_FourShells_Uneven_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Uneven_450_pc = (distance * x_FourShells_Uneven_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Uneven_450 = Stellar_Data_FourShells_Uneven_450['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Uneven_450 = stellar_FourShells_Uneven_450 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Even_450 = Table.read('Model_FourShells_EvenMass_JCMT_SCUBA2_450_RadialProfile.csv', format='ascii.csv')
x_FourShells_Even_450 = Stellar_Data_FourShells_Even_450['x_axis(arcsec)']
#x_FourShells_Even_450_cm = (distance * x_FourShells_Even_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Even_450_pc = (distance * x_FourShells_Even_450) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Even_450 = Stellar_Data_FourShells_Even_450['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Even_450 = stellar_FourShells_Even_450 * 1000 #Converting from Jy to mJy

############# Model Profiles for 850 micron (JCMT SCUBA 2 850) ##################################################

Stellar_Data_InnerShell_850 = Table.read('Model_InnerShell_Only_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_InnerShell_850 = Stellar_Data_InnerShell_850['x_axis(arcsec)']
#x_InnerShell_850_cm = (distance * x_InnerShell_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_InnerShell_850_pc = (distance * x_InnerShell_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_InnerShell_850 = Stellar_Data_InnerShell_850['Stellar_Profile(Jy/arcsec^2)']
stellar_InnerShell_850 = stellar_InnerShell_850 * 1000 #Converting from Jy to mJy

Stellar_Data_ShellTwo_850 = Table.read('Model_ShellTwo_Only_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_ShellTwo_850 = Stellar_Data_ShellTwo_850['x_axis(arcsec)']
#x_ShellTwo_850_cm = (distance * x_ShellTwo_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_ShellTwo_850_pc = (distance * x_ShellTwo_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_ShellTwo_850 = Stellar_Data_ShellTwo_850['Stellar_Profile(Jy/arcsec^2)']
stellar_ShellTwo_850 = stellar_ShellTwo_850 * 1000 #Converting from Jy to mJy

Stellar_Data_ShellThree_850 = Table.read('Model_ShellThree_Only_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_ShellThree_850 = Stellar_Data_ShellThree_850['x_axis(arcsec)']
#x_ShellThree_850_cm = (distance * x_ShellThree_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_ShellThree_850_pc = (distance * x_ShellThree_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_ShellThree_850 = Stellar_Data_ShellThree_850['Stellar_Profile(Jy/arcsec^2)']
stellar_ShellThree_850 = stellar_ShellThree_850 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_M2010_850 = Table.read('Model_OuterShell_Only_M2010_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_OuterShell_M2010_850 = Stellar_Data_OuterShell_M2010_850['x_axis(arcsec)']
#x_OuterShell_M2010_850_cm = (distance * x_OuterShell_M2010_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_M2010_850_pc = (distance * x_OuterShell_M2010_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_M2010_850 = Stellar_Data_OuterShell_M2010_850['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_M2010_850 = stellar_OuterShell_M2010_850 * 1000 #Converting from Jy to mJy

Stellar_Data_OuterShell_K2010_850 = Table.read('Model_OuterShell_Only_K2010_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_OuterShell_K2010_850 = Stellar_Data_OuterShell_K2010_850['x_axis(arcsec)']
#x_OuterShell_K2010_850_cm = (distance * x_OuterShell_K2010_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_OuterShell_K2010_850_pc = (distance * x_OuterShell_K2010_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_OuterShell_K2010_850 = Stellar_Data_OuterShell_K2010_850['Stellar_Profile(Jy/arcsec^2)']
stellar_OuterShell_K2010_850 = stellar_OuterShell_K2010_850 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Uneven_850 = Table.read('Model_FourShells_UnEvenMass_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_FourShells_Uneven_850 = Stellar_Data_FourShells_Uneven_850['x_axis(arcsec)']
#x_FourShells_Uneven_850_cm = (distance * x_FourShells_Uneven_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Uneven_850_pc = (distance * x_FourShells_Uneven_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Uneven_850 = Stellar_Data_FourShells_Uneven_850['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Uneven_850 = stellar_FourShells_Uneven_850 * 1000 #Converting from Jy to mJy

Stellar_Data_FourShells_Even_850 = Table.read('Model_FourShells_EvenMass_JCMT_SCUBA2_850_RadialProfile.csv', format='ascii.csv')
x_FourShells_Even_850 = Stellar_Data_FourShells_Even_850['x_axis(arcsec)']
#x_FourShells_Even_850_cm = (distance * x_FourShells_Even_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
x_FourShells_Even_850_pc = (distance * x_FourShells_Even_850) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
stellar_FourShells_Even_850 = Stellar_Data_FourShells_Even_850['Stellar_Profile(Jy/arcsec^2)']
stellar_FourShells_Even_850 = stellar_FourShells_Even_850 * 1000 #Converting from Jy to mJy




##### Reduced ChiSqaured Generation #################

print("Wavelength and Model Type, Reduced ChiSqaured Result") #Prints a heading in the file given when running python. 


#################### Reduced RedChiSqaured for 70 micron ##############################
radial_limit_obs_70 = x_70 < 100
radial_limit_Mod_70 = x_InnerShell_70 < 100

#Model Inner Shell Only
InnerShell_70_RedChiSq = ( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_InnerShell_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) )
                                                                                                                                / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
print("Model_InnerShell_Only_70,", InnerShell_70_RedChiSq) 

#Model Shell Two Only
ShellTwo_70_RedChiSq = ( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_ShellTwo_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) )
                                                                                                                                / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
print("Model_ShellTwo_Only_70,", ShellTwo_70_RedChiSq) 

#Model Shell Three Only
ShellThree_70_RedChiSq = ( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_ShellThree_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) )
                                                                                                                                / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
print("Model_ShellThree_Only_70,", ShellThree_70_RedChiSq) 

#Model Outer Shell (M2010) - Shell size from Maercker+2010
OuterShell_M2010_70_RedChiSq = ( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_OuterShell_M2010_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) ) 
                                                                                                                                                / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
#print(np.shape(stellar_70[radial_limit_obs_70]))
print("Model_OuterShell_Only_M2010_70,", OuterShell_M2010_70_RedChiSq) 

#Model Outer Shell (K2010) - Shell size from Kerschbaum+2010
OuterShell_K2010_70_RedChiSq =( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_OuterShell_K2010_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) )
                                                                                                                                    / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
print("Model_OuterShell_Only_K2010_70,", OuterShell_K2010_70_RedChiSq) 

#Model Four Shells - Dust Mass Unevenly spread between four shells
FourShells_Uneven_70_RedChiSq = ( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_FourShells_Uneven_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) )
                                                                                                                                    / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
print("Model_FourShells_UnEvenMass_70,", FourShells_Uneven_70_RedChiSq) 

#Model Four Shells + Cont Emission - Dust Mass Evenly spread between four shells
FourShells_Even_70_RedChiSq =( ( np.sum ( ( (stellar_70[radial_limit_obs_70] - stellar_FourShells_Even_70[radial_limit_Mod_70]) / stellar_unc_70[radial_limit_obs_70]) **2 ) )
                                                                                                                                        / ( np.shape(stellar_70[radial_limit_obs_70])[0] ) )
print("Model_FourShells_EvenMass_70,", FourShells_Even_70_RedChiSq) 



############# Reduced RedChiSqaured for 160 micron #############################
radial_limit_obs_160 = x_160 < 100
radial_limit_Mod_160 = x_InnerShell_160 < 100

#Model Inner Shell Only
InnerShell_160_RedChiSq = ( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_InnerShell_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) )
                                                                                                                                / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
print("Model_InnerShell_Only_160,", InnerShell_160_RedChiSq) 

#Model Shell Two Only
ShellTwo_160_RedChiSq = ( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_ShellTwo_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) )
                                                                                                                                / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
print("Model_ShellTwo_Only_160,", ShellTwo_160_RedChiSq) 

#Model Shell Three Only
ShellThree_160_RedChiSq = ( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_ShellThree_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) )
                                                                                                                                / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
print("Model_ShellThree_Only_160,", ShellThree_160_RedChiSq) 

#Model Outer Shell (M2010) - Shell size from Maercker+2010
OuterShell_M2010_160_RedChiSq = ( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_OuterShell_M2010_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) ) 
                                                                                                                                                / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
#print(np.shape(stellar_160[radial_limit_obs_160]))
print("Model_OuterShell_Only_M2010_160,", OuterShell_M2010_160_RedChiSq) 

#Model Outer Shell (K2010) - Shell size from Kerschbaum+2010
OuterShell_K2010_160_RedChiSq =( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_OuterShell_K2010_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) )
                                                                                                                                    / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
print("Model_OuterShell_Only_K2010_160,", OuterShell_K2010_160_RedChiSq) 

#Model Four Shells - Dust Mass Unevenly spread between four shells
FourShells_Uneven_160_RedChiSq = ( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_FourShells_Uneven_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) )
                                                                                                                                    / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
print("Model_FourShells_UnEvenMass_160,", FourShells_Uneven_160_RedChiSq) 

#Model Four Shells - Dust Mass Evenly spread between four shells
FourShells_Even_160_RedChiSq =( ( np.sum ( ( (stellar_160[radial_limit_obs_160] - stellar_FourShells_Even_160[radial_limit_Mod_160]) / stellar_unc_160[radial_limit_obs_160]) **2 ) )
                                                                                                                                        / ( np.shape(stellar_160[radial_limit_obs_160])[0] ) )
print("Model_FourShells_EvenMass_160,", FourShells_Even_160_RedChiSq) 



############# Reduced RedChiSqaured for 450 micron ######################
radial_limit_obs_450 = x_450 < 60
radial_limit_Mod_450 = x_InnerShell_450 < 60

#Model Inner Shell Only
InnerShell_450_RedChiSq = ( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_InnerShell_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) )
                                                                                                                                / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) )
print("Model_InnerShell_Only_450,", InnerShell_450_RedChiSq) 

#Model Shell Two Only
ShellTwo_450_RedChiSq = ( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_ShellTwo_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) )
                                                                                                                                / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) )
print("Model_ShellTwo_Only_450,", ShellTwo_450_RedChiSq) 

#Model Shell Three Only
ShellThree_450_RedChiSq = ( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_ShellThree_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) )
                                                                                                                                / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) )
print("Model_ShellThree_Only_450,", ShellThree_450_RedChiSq) 

#Model Outer Shell (M2010) - Shell size from Maercker+2010
OuterShell_M2010_450_RedChiSq = ( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_OuterShell_M2010_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) ) 
                                                                                                                                                / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) )
#print(np.shape(stellar_450[radial_limit_obs_450]))
print("Model_OuterShell_Only_M2010_450,", OuterShell_M2010_450_RedChiSq) 

#Model Outer Shell (K2010) - Shell size from Kerschbaum+2010
OuterShell_K2010_450_RedChiSq =( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_OuterShell_K2010_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) )
                                                                                                                                    / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) )
print("Model_OuterShell_Only_K2010_450,", OuterShell_K2010_450_RedChiSq) 

#Model Four Shells - Dust Mass Unevenly spread between four shells
FourShells_Uneven_450_RedChiSq = ( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_FourShells_Uneven_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) )
                                                                                                                                    / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) )
print("Model_FourShells_UnEvenMass_450,", FourShells_Uneven_450_RedChiSq) 

#Model Four Shells Dust Mass Evenly spread between four shells
FourShells_Even_450_RedChiSq =( ( np.sum ( ( (stellar_450[radial_limit_obs_450] - stellar_FourShells_Even_450[radial_limit_Mod_450]) / stellar_unc_450[radial_limit_obs_450]) **2 ) )
                                                                                                                                        / ( np.shape(stellar_450[radial_limit_obs_450])[0] ) )
print("Model_FourShells_EvenMass_450,", FourShells_Even_450_RedChiSq) 





######################## Reduced RedChiSqaured for 850 micron ############################
radial_limit_obs_850 = x_850 < 100
radial_limit_Mod_850 = x_InnerShell_850 < 100

#Model Inner Shell Only
InnerShell_850_RedChiSq = ( ( np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_InnerShell_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 ) )
                                                                                                                                / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
print("Model_InnerShell_Only_850,", InnerShell_850_RedChiSq) 

#Model Shell Two Only
ShellTwo_850_RedChiSq = ( ( np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_ShellTwo_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 ) )
                                                                                                                                / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
print("Model_ShellTwo_Only_850,", ShellTwo_850_RedChiSq) 

#Model Shell Three Only
ShellThree_850_RedChiSq = ( ( np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_ShellThree_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 ) )
                                                                                                                                / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
print("Model_ShellThree_Only_850,", ShellThree_850_RedChiSq) 

#Model Outer Shell (M2010) - Shell size from Maercker+2010
OuterShell_M2010_850_RedChiSq = ( ( np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_OuterShell_M2010_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 ) ) 
                                                                                                                                                / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
#print(np.shape(stellar_850[radial_limit_obs_850]))
print("Model_OuterShell_Only_M2010_850,", OuterShell_M2010_850_RedChiSq) 

#Model Outer Shell (K2010) - Shell size from Kerschbaum+2010
OuterShell_K2010_850_RedChiSq =( ( np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_OuterShell_K2010_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 ) )
                                                                                                                                    / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
print("Model_OuterShell_Only_K2010_850,", OuterShell_K2010_850_RedChiSq) 

#Model Four Shells - Dust Mass Unevenly spread between four shells
FourShells_Uneven_850_RedChiSq = ( ( np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_FourShells_Uneven_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 ) )
                                                                                                                                    / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
print("Model_FourShells_UnEvenMass_850,", FourShells_Uneven_850_RedChiSq) 

#Model Four Shells Dust Mass Evenly spread between four shells
FourShells_Even_850_RedChiSq =( ( np.sum ( ( (stellar_850[radial_limit_obs_850] - stellar_FourShells_Even_850[radial_limit_Mod_850]) / stellar_unc_850[radial_limit_obs_850]) **2 ) )
                                                                                                                                        / ( np.shape(stellar_850[radial_limit_obs_850])[0] ) )
print("Model_FourShells_EvenMass_850,", FourShells_Even_850_RedChiSq) 
















































