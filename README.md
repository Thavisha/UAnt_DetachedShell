# U Ant Detached Dust Shell in the sub-mm
Reproducing "The Nearby Evolved Stars Survey: I. JCMT/SCUBA-2Sub-millimetre detection of the detached shell of UAntliae" paper by Thavisha E. Dharmawardena, et al., subm. MNRAS. 

A full description of the methods used is provided in Secs. 2 and 3 in the paper.  

Generated using python3.6 (suitable for python2.7 and later versions)

Steps followed:
1) Data Download and Reduction
2) Radial profile generation and plotting
3) Aperture Photometry 
4) Generating images of the SCUBA-2 and HARP observations
5) Radial point-to-point SED fitting and determining Temperature, Density and Beta profiles using EMCEE
6) Calculating mass-loss from the Density profiles
7) Raditative transfer modelling using Hyperion


For testing the differences between the CO-subtracted and Non-CO-subtracted SCUBA-2 850micron observations replace the CO-subtracted .fits file with the non-CO-subtracted .fits file and use the same code.



