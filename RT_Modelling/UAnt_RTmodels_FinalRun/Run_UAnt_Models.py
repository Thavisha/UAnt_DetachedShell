import numpy as np
from hyperion.model import AnalyticalYSOModel
from hyperion.util.constants import rsun, lsun, au, msun, year, c, sigma #sigma = notation for stefan-boltzman constant

from Model_OuterShell_ContinousEmission_M2010 import OuterShell_plusCont_M2010
from Model_ContinousEmission_Only import Cont_Only
from Model_InnerShell_ContinousEmission import InnerShell_plusCont
from Model_FourShells_ContinousEmission_EvenMass import FourShells_plusCont_EvenMass
from Model_FourShells_ContinousEmission_UnevenMass import FourShells_plusCont_UnevenMass
from Model_OuterShell_ContinousEmission_K2010 import OuterShell_plusCont_K2010

"""
Single file to call and run all different models set up for U Ant.

(It is fine to leave the continous emission going from 1.5R* to 80" because it only contributes ~1.5% of the total mass if the detached shells are present. Therefore if the detached shells are present on top of the cont. emission their mass will dominate by a lot and the small 1.5% mass contribution underneath will not matter. We do not need to break the continious emission to pre and post therm.pulse.)
 
"""


#Input Source Parameters
uant_distance = 268 #pc
stellar_spectrum_name = 'Aringer+2009_StellarSpectrum.ascii' #File Cols = wav-Angstroms, fcont-NotNeeded-Ignore, nufnu-Convert to fnu 
dust_file_name = 'Cstar_DustProperties.hdf5'   #'hyperion-dust-0.1.0/dust_files/kmh94_3.1_full.hdf5'
total_shell_mass = 2e-5 * msun #From my SED fit (point-point) Integrated under the Sigma (density) profile - Fitted after CO subtraction from SCUBA-2 850 
expansion_velocity = (4.5 * 10e5) * year #CO expsansion vel averaged from Bands6&3 from Kerschbaum+2017 converted from 4.5 km/s to cm/year for easy conversion to derive total constant outflow mass later on. 
present_day_mlr = (3.8e-12) * 0.45 #From GRAMS fit. #* 0.45 - GRAMS measured the the present day dust mlr for an exp.vel of 10 km/s. As our exp.vel is 4.5 km/s we must scale the GRAMS val. to match that. 
stellar_luminoisity = 6920 * lsun #From GRAMS fit.
stellar_temperature = 2600 #Units=Kelvin #From GRAMS fit.
stellar_radius = 424 * rsun #From GRAMS fit.
continuous_envelope_inner_rad = 1.5 * stellar_radius #From GRAMS fit.
continuous_envelope_outer_rad = 80 #Units=arcsec. Outer radius of continuous envelope - estimating the total emission will be ~ this radius taking in to account pre therm.pulse mass loss. From my SED fit (point-point) Integrated under the Sigma (density) profile
 

Photon_Number = 1e9  #Input Number of Photons. #Number of photons to run the RT models with (Short run - 1e7 // Long run - 1e9)
CPU_Number = 30 #Number of CPUs to use when running hyperion (8 for my laptop. Up to 32 for Baldrick - I use 30)


# Read in stellar spectrum for Aringer+2009_StellarSpectrum.ascii
wav, fcont, nufnu = np.loadtxt(stellar_spectrum_name, unpack=True) #wav, fcont, nufnu. Hyperion needs nu and fnu so we convert out file values to this
nu = c / (wav * 1.e-8) #(wav * 1.e-8) to convert from angstrom to cm
fnu = nufnu / nu

extraploation_nu = 3e11 #Extrapoplate up to 3e8 Hz (3e8=100 cm) (3e11=1mm)to have the stellar spectrum at long wavelengths 
extrapolation_function = ((extraploation_nu/nu[-1]) ** 2 ) * fnu[-1]
fnu = np.append(fnu, extrapolation_function)
nu = np.append(nu, extraploation_nu)
#print(nu, fnu)




#Model 1 = Outer Shell + Continous Emission Model (Using Radii from Maercker et al., 2010 - Scattered light Data - Outermost shell - shell 4)
Input_OuterShell_ContEm_M2010, Output_OuterShell_ContEm_M2010 = OuterShell_plusCont_M2010(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature, stellar_radius, total_shell_mass, uant_distance, present_day_mlr, expansion_velocity, continuous_envelope_inner_rad, continuous_envelope_outer_rad, Photon_Number, CPU_Number)

#Model 2 = Continous Emission Only Model
Input_ContEm_Only, Output_ContEm_Only = Cont_Only(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature, stellar_radius, total_shell_mass, uant_distance, present_day_mlr, expansion_velocity, continuous_envelope_inner_rad,  continuous_envelope_outer_rad, Photon_Number, CPU_Number)

#Model 3 = Inner Shell + Continous Emission Model
Input_InnerShell_ContEm, Output_InnerShell_ContEm = InnerShell_plusCont(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature, stellar_radius, total_shell_mass, uant_distance, present_day_mlr, expansion_velocity, continuous_envelope_inner_rad, continuous_envelope_outer_rad, Photon_Number, CPU_Number)

#Model 4 = Four Shells + Continious Emission Model - Shell Mass distributed Evenly
Input_FourShells_ContEm_EvenMass, Output_FourShells_ContEm_EvenMass = FourShells_plusCont_EvenMass(dust_file_name, fnu, nu, stellar_luminoisity,  stellar_radius, stellar_temperature, total_shell_mass, uant_distance, present_day_mlr, expansion_velocity, continuous_envelope_inner_rad, continuous_envelope_outer_rad, Photon_Number, CPU_Number)

#Model 5 = Four Shells + Continious Emission Model - Shell Mass distributed Unevenly in the four shells 50% in shell 4, 5% in shell 3, 23% in shell 2, 23% in shell 1 
Input_FourShells_ContEm_UnevenMass, Output_FourShells_ContEm_UnevenMass = FourShells_plusCont_UnevenMass(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature, stellar_radius, total_shell_mass, uant_distance, present_day_mlr, expansion_velocity, continuous_envelope_inner_rad, continuous_envelope_outer_rad, Photon_Number, CPU_Number)

#Model 6 = Outer Shell + Continous Emission Model (Using Radii from Kerschbaum et al., 2010 - PACS data - Outermost shell - shell 3)
Input_OuterShell_ContEm_K2010, Output_OuterShell_ContEm_K2010 = OuterShell_plusCont_K2010(dust_file_name, fnu, nu, stellar_luminoisity, stellar_temperature, stellar_radius, total_shell_mass, uant_distance, present_day_mlr, expansion_velocity, continuous_envelope_inner_rad, continuous_envelope_outer_rad, Photon_Number, CPU_Number)


