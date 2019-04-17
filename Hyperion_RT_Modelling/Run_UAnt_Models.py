import numpy as np
from hyperion.model import AnalyticalYSOModel
from hyperion.util.constants import rsun, lsun, au, msun, year, c, sigma #sigma = notation for stefan-boltzman constant


from Model_InnerShell_Only import InnerShell_Only
from Model_ShellTwo_Only import ShellTwo_Only
from Model_ShellThree_Only import ShellThree_Only
from Model_OuterShell_Only_M2010 import OuterShell_Only_M2010
from Model_OuterShell_Only_K2010 import OuterShell_Only_K2010
from Model_FourShells_UnEvenMass import FourShells_UnEvenMass
from Model_FourShells_EvenMass import FourShells_EvenMass





"""
Single file to call and run all different models set up for U Ant.

It's fine to ignore the present day MLR or the continuous pre and post TP mass loss since it's so tiny compared to the 
"""


#Input Source Parameters
uant_distance = 268 #pc
stellar_spectrum_name = 'Aringer+2009_StellarSpectrum.ascii' #File Cols = wav-Angstroms, fcont-NotNeeded-Ignore, nufnu-Convert to fnu 
dust_file_name = 'Cstar_DustProperties.hdf5'   #'hyperion-dust-0.1.0/dust_files/kmh94_3.1_full.hdf5'
total_shell_mass = 2e-5 * msun #From my SED fit (point-point) Integrated under the Sigma (density) profile - Fitted after CO subtraction from SCUBA-2 850 
expansion_velocity = (4.5 * 10e5) * year #CO expansion vel averaged from Bands6&3 from Kerschbaum+2017 converted from 4.5 km/s to cm/year for easy conversion to derive total constant outflow mass later on.  
stellar_luminoisity = 6920 * lsun #L of U Ant fit.
stellar_temperature = 2600 #Units=Kelvin #From COMARCS fit.
stellar_radius = 424 * rsun #From COMARCS fit.
MaxCSE_outer_rad = 80 #Units=arcsec. Outer radius of CSE - From my SED fit (point-point) Integrated under the Sigma (density) profile
 

Photon_Number = 1e9  #Input Number of Photons. #Number of photons to run the RT models with (Short run - 1e6 // Long run - 1e9)
CPU_Number = 30 #Number of CPUs to use when running Hyperion (8 for my laptop. Up to 32 for Baldrick - I use 30)


# Read in stellar spectrum for Aringer+2009_StellarSpectrum.ascii
wav, fcont, nufnu = np.loadtxt(stellar_spectrum_name, unpack=True) #wav, fcont, nufnu. Hyperion needs nu and fnu so we convert out file values to this
nu = c / (wav * 1.e-8) #(wav * 1.e-8) to convert from angstrom to cm
fnu = nufnu / nu #In reality Aringer Model gives nuLnu but the difference between that and nuFnu is scaling factor which Hyperion automatically account for by scaling up to the given stellar luminosity so we don't need worry about converting from one to the other

extraploation_nu = 3e11 #Extrapoplate up to 3e8 Hz (3e8=100 cm) (3e11=1mm)to have the stellar spectrum at long wavelengths 
extrapolation_function = ((extraploation_nu/nu[-1]) ** 2 ) * fnu[-1]
fnu = np.append(fnu, extrapolation_function)
nu = np.append(nu, extraploation_nu)
#print(nu, fnu)




#Model 1: Inner shell only  (All shell dust mass put into shell 1 defined by G2001&2003)
Input_InnerShell_Only, Output_InnerShell_Only = InnerShell_Only(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature,  stellar_radius, total_shell_mass, uant_distance, expansion_velocity, MaxCSE_outer_rad, Photon_Number, CPU_Number)

#Model 2: Shell Two Only (All shell dust mass put into shell 2 defined by G2001&2003) 
Input_ShellTwo_Only, Output_ShellTwo_Only = ShellTwo_Only(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature,  stellar_radius, total_shell_mass, uant_distance, expansion_velocity, MaxCSE_outer_rad, Photon_Number, CPU_Number)

#Model 3: Shell Three Only (All shell dust mass put into shell 2 defined by G2001&2003 AND M2010) 
Input_ShellThree_Only, Output_ShellThree_Only = ShellThree_Only(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature,  stellar_radius, total_shell_mass, uant_distance, expansion_velocity, MaxCSE_outer_rad, Photon_Number, CPU_Number)

#Model 4: Shell 4 Only with Radius of shell as given in M2010 - All dust mass in shell 4
Input_OuterShell_Only_M2010, Output_OuterShell_Only_M2010 = OuterShell_Only_M2010(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature,  stellar_radius, total_shell_mass, uant_distance, expansion_velocity, MaxCSE_outer_rad, Photon_Number, CPU_Number)

#Model 5: Shell 4 Only with Radius of shell as given in K2010 - All dust mass in shell 4
Input_OuterShell_Only_K2010, Output_OuterShell_Only_K2010 = OuterShell_Only_K2010(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature,  stellar_radius, total_shell_mass, uant_distance, expansion_velocity, MaxCSE_outer_rad, Photon_Number, CPU_Number)

#Model 6: Four shells Model - Shell mass distributed unevenly based on literature values
Input_FourShells_UnevenMass, Output_FourShells_UnevenMass = FourShells_UnEvenMass(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature, stellar_radius, total_shell_mass, uant_distance, expansion_velocity,  MaxCSE_outer_rad, Photon_Number, CPU_Number)

#Model 7: Four Shells Model - Shell Mass distributed Evenly
Input_FourShells_EvenMass, Output_FourShells_EvenMass = FourShells_EvenMass(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature,  stellar_radius, total_shell_mass, uant_distance, expansion_velocity, MaxCSE_outer_rad, Photon_Number, CPU_Number)
