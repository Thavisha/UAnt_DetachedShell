import matplotlib.pyplot as plt

from hyperion.model import ModelOutput
from hyperion.util.constants import pc

m = ModelOutput('AGB_Model_Example_6.rtout')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Extract all SEDs
sed = m.get_sed(inclination='all', aperture=-1, distance=268 * pc, units='Jy')
#help(sed) #Needed to check the SED file headers. 



# Plot SED for each inclination
ax.errorbar(sed.wav, sed.val, yerr=sed.unc, color='black')

ax.set_xscale('log') #setting x axis to log scale
ax.set_yscale('log') #setting y axis to log scale
ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel(r'$F_\nu$ [Jy]')
ax.set_xlim(0.1, 1500.)
#ax.set_ylim(2.e-16, 2.e-9)
fig.savefig('AGB_Model_Example_5.png')

plt.show()
