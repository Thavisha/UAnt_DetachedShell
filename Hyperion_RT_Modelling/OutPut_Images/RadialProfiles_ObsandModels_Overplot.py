import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
#import pylatex
import matplotlib.ticker as mtick
import warnings

"""
Plotting RT Models and Observed Radial profiles together for Comparison

Code modified from Obsevred Radial profile plotting

#Heavily Edited for U Ant RT Models and Obs Plotting!! Don't use for anything else!!!!!!!!!!

"""

font = {'family' : 'normal',
        'size'   : 24,
	'weight' : 'medium'}



plt.rc('font', **font)

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





#################### Plotting ##################################################################


fig = plt.figure(figsize=(20, 20)) #(width,height)
gs = gridspec.GridSpec(4, 1) #(no.of.row, no.of.cols, height.ratio.of.rows)
#gs = gridspec.GridSpec(2, 2)
gs.update(hspace=0)
gs.update(wspace=0)
#length = [100, 60, 30, 60] #x axis arcseconds lengths for 70, 160, 450, 850 in this specific order!
#b = x[i] < length[i]	

#Adding common x and y labels
fig.text(0.52, 0.07,"Radius ($^{\prime\prime}$)", ha='center') #(distance.from.left.edge.of.image, distance.from.bottom.of.image)
#fig.text(0.52, 0.93, 'Radius ($\\times 10^{18}$cm)', ha='center') #cm units
fig.text(0.52, 0.91, 'Radius (pc)', ha='center') #pc units
fig.text(0.06, 0.5,"Normalised Surface Brightness", va='center', rotation='vertical')


#Plotting 70micron Radial + PSF Profile
ax1 = fig.add_subplot(gs[0,0]) #gs[row, col]
radial_limit_obs_70 = x_70 < 100

#We plot the normalised surface brightness profile (i.e: profile/profile.max) so then we need the normalised unc for the obs. prof as well.
frac_unc_70 =np.sqrt(
                ( ( stellar_unc_70[radial_limit_obs_70] / stellar_70[radial_limit_obs_70] ) ** 2 )
                    +( ( stellar_unc_70[np.argmax(stellar_70[radial_limit_obs_70])] / np.max(stellar_70[radial_limit_obs_70]) )**2 ) 
                    )
norm_stellar_70 = stellar_70[radial_limit_obs_70]/np.max(stellar_70[radial_limit_obs_70])
norm_stellar_unc_70 = frac_unc_70 * norm_stellar_70

ax1.errorbar(x_70[radial_limit_obs_70], norm_stellar_70, yerr=norm_stellar_unc_70, fmt='o', color='black', markersize=8, capsize=7, label='70$\\mu m$ Observed Profile')
radial_limit_Mod_70 = x_ContEm_Only_70 < 100
ax1.errorbar(x_ContEm_Only_70[radial_limit_Mod_70], stellar_ContEm_Only_70[radial_limit_Mod_70]/np.max(stellar_ContEm_Only_70[radial_limit_Mod_70]), fmt='-', color='violet', markersize=8, capsize=7, label='70$\\mu m$ Model_NoShells') 
ax1.errorbar(x_ContEm_Only_70[radial_limit_Mod_70], stellar_InnerShell_70[radial_limit_Mod_70]/np.max(stellar_InnerShell_70[radial_limit_Mod_70]), fmt='-', color='green', markersize=8, capsize=7, label='70$\\mu m$ Model_Inner')
ax1.errorbar(x_ContEm_Only_70[radial_limit_Mod_70], stellar_OuterShell_K2010_70[radial_limit_Mod_70]/np.max(stellar_OuterShell_K2010_70[radial_limit_Mod_70]), fmt='-', color='blue', markersize=8, capsize=7, label='70$\\mu m$ Model_Outer_K2010')
ax1.errorbar(x_ContEm_Only_70[radial_limit_Mod_70], stellar_OuterShell_M2010_70[radial_limit_Mod_70]/np.max(stellar_OuterShell_M2010_70[radial_limit_Mod_70]), fmt='-', color='red', markersize=8, capsize=7, label='70$\\mu m$ Model_Outer_M2010')
ax1.errorbar(x_ContEm_Only_70[radial_limit_Mod_70], stellar_FourShells_Even_70[radial_limit_Mod_70]/np.max(stellar_FourShells_Even_70[radial_limit_Mod_70]), fmt='-', color='purple', markersize=8, capsize=7, label='70$\\mu m$ Model_FourShells_EvenMass')
ax1.errorbar(x_ContEm_Only_70[radial_limit_Mod_70], stellar_FourShells_Uneven_70[radial_limit_Mod_70]/np.max(stellar_FourShells_Uneven_70[radial_limit_Mod_70]), fmt='-', color='orange', markersize=8, capsize=7, label='70$\\mu m$ Model_FourShells_UnevenMass')

nbins = len(ax1.get_yticklabels())
ax1.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both'))
ax1Ticks = ax1.get_xticks() 
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_yscale('log') #setting y axis to log scale
#ax1.set_ylim(0.02, 1.2)
#ax1.legend(fontsize=11)
ax1.set_xlim(xmin=0) #starting x axis from zero

#Adding text box with wavelength label for PACS 70
textstr = 'Hershcel/PACS 70 $\mu$m'
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax1.text(0.98, 0.9, textstr, transform=ax1.transAxes, fontsize=18, verticalalignment='top', horizontalalignment='right', bbox=props) #(xPos, yPos, text, ....)





#Plotting 160micron Radial + PSF Profile
ax2 = fig.add_subplot(gs[1,0], sharex=ax1) #gs[row, col]
radial_limit_obs_160 = x_160 < 100

frac_unc_160 =np.sqrt(
                ( ( stellar_unc_160[radial_limit_obs_160] / stellar_160[radial_limit_obs_160] ) ** 2 )
                    +( ( stellar_unc_160[np.argmax(stellar_160[radial_limit_obs_160])] / np.max(stellar_160[radial_limit_obs_160]) )**2 ) 
                    )
norm_stellar_160 = stellar_160[radial_limit_obs_160]/np.max(stellar_160[radial_limit_obs_160])
norm_stellar_unc_160 = frac_unc_160 * norm_stellar_160

ax2.errorbar(x_160[radial_limit_obs_160], norm_stellar_160, yerr=norm_stellar_unc_160, fmt='o', color='black', markersize=8, capsize=7, label='160$\\mu m$ Observed Profile')

radial_limit_Mod_160 = x_ContEm_Only_160 < 100
ax2.errorbar(x_ContEm_Only_160[radial_limit_Mod_160], stellar_ContEm_Only_160[radial_limit_Mod_160]/np.max(stellar_ContEm_Only_160[radial_limit_Mod_160]), fmt='-', color='violet', markersize=8, capsize=7, label='160$\\mu m$ Model_NoShells') 
ax2.errorbar(x_ContEm_Only_160[radial_limit_Mod_160], stellar_InnerShell_160[radial_limit_Mod_160]/np.max(stellar_InnerShell_160[radial_limit_Mod_160]), fmt='-', color='green', markersize=8, capsize=7, label='160$\\mu m$ Model_Inner')
ax2.errorbar(x_ContEm_Only_160[radial_limit_Mod_160], stellar_OuterShell_K2010_160[radial_limit_Mod_160]/np.max(stellar_OuterShell_K2010_160[radial_limit_Mod_160]), fmt='-', color='blue', markersize=8, capsize=7, label='160$\\mu m$ Model_Outer_K2010')
ax2.errorbar(x_ContEm_Only_160[radial_limit_Mod_160], stellar_OuterShell_M2010_160[radial_limit_Mod_160]/np.max(stellar_OuterShell_M2010_160[radial_limit_Mod_160]), fmt='-', color='red', markersize=8, capsize=7, label='160$\\mu m$ Model_Outer_M2010')
ax2.errorbar(x_ContEm_Only_160[radial_limit_Mod_160], stellar_FourShells_Even_160[radial_limit_Mod_160]/np.max(stellar_FourShells_Even_160[radial_limit_Mod_160]), fmt='-', color='purple', markersize=8, capsize=7, label='160$\\mu m$ Model_FourShells_EvenMass')
ax2.errorbar(x_ContEm_Only_160[radial_limit_Mod_160], stellar_FourShells_Uneven_160[radial_limit_Mod_160]/np.max(stellar_FourShells_Uneven_160[radial_limit_Mod_160]), fmt='-', color='orange', markersize=8, capsize=7, label='160$\\mu m$ Model_FourShells_UnevenMass') 

nbins = len(ax2.get_yticklabels())
ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both'))
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_yscale('log') #setting y axis to log scale
#ax2.legend(fontsize=11, loc=1)
ax2.set_xlim(xmin=0) #starting x axis from zero

#Adding text box with wavelength label for PACS 160
textstr = 'Hershcel/PACS 160 $\mu$m'
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax2.text(0.98, 0.9, textstr, transform=ax2.transAxes, fontsize=18, verticalalignment='top', horizontalalignment='right', bbox=props) #(xPos, yPos, text, ....)




#Plotting 450micron Radial + PSF Profile
ax3 = fig.add_subplot(gs[2,0], sharex=ax1) #gs[row, col]
radial_limit_obs_450 = x_450 < 60

frac_unc_450 =np.sqrt(
                ( ( stellar_unc_450[radial_limit_obs_450] / stellar_450[radial_limit_obs_450] ) ** 2 )
                    +( ( stellar_unc_450[np.argmax(stellar_450[radial_limit_obs_450])] / np.max(stellar_450[radial_limit_obs_450]) )**2 ) 
                    )
norm_stellar_450 = stellar_450[radial_limit_obs_450]/np.max(stellar_450[radial_limit_obs_450])
norm_stellar_unc_450 = frac_unc_450 * norm_stellar_450

ax3.errorbar(x_450[radial_limit_obs_450], norm_stellar_450, yerr=norm_stellar_unc_450, fmt='o', color='black', markersize=8, capsize=7, label='450$\\mu m$ Observed Profile')

radial_limit_Mod_450 = x_ContEm_Only_450 < 60
ax3.errorbar(x_ContEm_Only_450[radial_limit_Mod_450], stellar_ContEm_Only_450[radial_limit_Mod_450]/np.max(stellar_ContEm_Only_450[radial_limit_Mod_450]), fmt='-', color='violet', markersize=8, capsize=7, label='450$\\mu m$ Model_NoShells') 
ax3.errorbar(x_ContEm_Only_450[radial_limit_Mod_450], stellar_InnerShell_450[radial_limit_Mod_450]/np.max(stellar_InnerShell_450[radial_limit_Mod_450]), fmt='-', color='green', markersize=8, capsize=7, label='450$\\mu m$ Model_Inner')
ax3.errorbar(x_ContEm_Only_450[radial_limit_Mod_450], stellar_OuterShell_K2010_450[radial_limit_Mod_450]/np.max(stellar_OuterShell_K2010_450[radial_limit_Mod_450]), fmt='-', color='blue', markersize=8, capsize=7, label='450$\\mu m$ Model_Outer_K2010')
ax3.errorbar(x_ContEm_Only_450[radial_limit_Mod_450], stellar_OuterShell_M2010_450[radial_limit_Mod_450]/np.max(stellar_OuterShell_M2010_450[radial_limit_Mod_450]), fmt='-', color='red', markersize=8, capsize=7, label='450$\\mu m$ Model_Outer_M2010')
ax3.errorbar(x_ContEm_Only_450[radial_limit_Mod_450], stellar_FourShells_Even_450[radial_limit_Mod_450]/np.max(stellar_FourShells_Even_450[radial_limit_Mod_450]), fmt='-', color='purple', markersize=8, capsize=7, label='450$\\mu m$ Model_FourShells_EvenMass')
ax3.errorbar(x_ContEm_Only_450[radial_limit_Mod_450], stellar_FourShells_Uneven_450[radial_limit_Mod_450]/np.max(stellar_FourShells_Uneven_450[radial_limit_Mod_450]), fmt='-', color='orange', markersize=8, capsize=7, label='450$\\mu m$ Model_FourShells_UnevenMass') 

nbins = len(ax3.get_yticklabels())
ax3.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both'))
plt.setp(ax3.get_xticklabels(), visible=False)
ax3.set_yscale('log') #setting y axis to log scale
ax3.legend(fontsize=15, loc=4)
ax3.set_xlim(xmin=0) #starting x axis from zero

#Adding text box with wavelength label for SCUBA2 450
textstr = 'JCMT/SCUBA-2 450 $\mu$m'
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax3.text(0.98, 0.9, textstr, transform=ax3.transAxes, fontsize=18, verticalalignment='top', horizontalalignment='right', bbox=props) #(xPos, yPos, text, ....)



#Plotting 850micron Radial + PSF Profile
ax4 = fig.add_subplot(gs[3,0], sharex=ax1) #gs[row, col]
radial_limit_obs_850 = x_850 < 100

frac_unc_850 =np.sqrt(
                ( ( stellar_unc_850[radial_limit_obs_850] / stellar_850[radial_limit_obs_850] ) ** 2 )
                    +( ( stellar_unc_850[np.argmax(stellar_850[radial_limit_obs_850])] / np.max(stellar_850[radial_limit_obs_850]) )**2 ) 
                    )
norm_stellar_850 = stellar_850[radial_limit_obs_850]/np.max(stellar_850[radial_limit_obs_850])
norm_stellar_unc_850 = frac_unc_850 * norm_stellar_850

ax4.errorbar(x_850[radial_limit_obs_850], norm_stellar_850, yerr=norm_stellar_unc_850, fmt='o', color='black', markersize=8, capsize=7, label='850$\\mu m$ Observed Profile')


radial_limit_Mod_850 = x_ContEm_Only_850 < 100
ax4.errorbar(x_ContEm_Only_850[radial_limit_Mod_850], stellar_ContEm_Only_850[radial_limit_Mod_850]/np.max(stellar_ContEm_Only_850[radial_limit_Mod_850]), fmt='-', color='violet', markersize=8, capsize=7, label='850$\\mu m$ Model_NoShells') 
ax4.errorbar(x_ContEm_Only_850[radial_limit_Mod_850], stellar_InnerShell_850[radial_limit_Mod_850]/np.max(stellar_InnerShell_850[radial_limit_Mod_850]), fmt='-', color='green', markersize=8, capsize=7, label='850$\\mu m$ Model_Inner')
ax4.errorbar(x_ContEm_Only_850[radial_limit_Mod_850], stellar_OuterShell_K2010_850[radial_limit_Mod_850]/np.max(stellar_OuterShell_K2010_850[radial_limit_Mod_850]), fmt='-', color='blue', markersize=8, capsize=7, label='850$\\mu m$ Model_Outer_K2010')
ax4.errorbar(x_ContEm_Only_850[radial_limit_Mod_850], stellar_OuterShell_M2010_850[radial_limit_Mod_850]/np.max(stellar_OuterShell_M2010_850[radial_limit_Mod_850]), fmt='-', color='red', markersize=8, capsize=7, label='850$\\mu m$ Model_Outer_M2010')
ax4.errorbar(x_ContEm_Only_850[radial_limit_Mod_850], stellar_FourShells_Even_850[radial_limit_Mod_850]/np.max(stellar_FourShells_Even_850[radial_limit_Mod_850]), fmt='-', color='purple', markersize=8, capsize=7, label='850$\\mu m$ Model_FourShells_EvenMass')
ax4.errorbar(x_ContEm_Only_850[radial_limit_Mod_850], stellar_FourShells_Uneven_850[radial_limit_Mod_850]/np.max(stellar_FourShells_Uneven_850[radial_limit_Mod_850]), fmt='-', color='orange', markersize=8, capsize=7, label='850$\\mu m$ Model_FourShells_UnevenMass') 

nbins = len(ax4.get_yticklabels())
ax4.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
ax4.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max val on x axis to prevent merging with res axis labels
ax4.set_yscale('log') #setting y axis to log scale
#ax4.legend(fontsize=11, loc=1)
ax4.set_xlim(xmin=0) #starting x axis from zero

#Adding text box with wavelength label for SCUBA2 850
textstr = 'JCMT/SCUBA-2 850 $\mu$m'
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax4.text(0.98, 0.9, textstr, transform=ax4.transAxes, fontsize=18, verticalalignment='top', horizontalalignment='right', bbox=props) #(xPos, yPos, text, ....)


#Adding second x axis on top of plot for the cm x axis scale using the cm scale for 70micron for both stellar and res profiles.
#radial_limit_obs_70[0:2]=True
ax5 = ax1.twiny()
ax5.errorbar(x_70_pc[radial_limit_obs_70], norm_stellar_70, yerr=norm_stellar_unc_70, fmt='o', color='none', mec='none') #pc units!
ax5.set_yscale('log') #setting y axis to log scale
ax5.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max val on x axis to prevent merging with res axis labels
ax5.set_xlim(xmin=0) #starting x axis from zero




#ax1.errorbar(x_70[radial_limit_obs_70], norm_stellar_70, yerr=norm_stellar_unc_70, fmt='o', color='black', markersize=8, capsize=7, label='70$\\mu m$ Observed Profile')


save_file = star + 'ObservedandModel_RadialProfile_withLegend.png'
#plt.savefig(save_file, bbox_inches='tight')
plt.show()




















