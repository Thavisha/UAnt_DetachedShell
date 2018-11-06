import numpy as np
from scipy.ndimage import zoom

"""
Rebinning Image to match scaling of observed images get it out at Jy/arcsec**2

"""



def Image_Rebin(BeamConv_Image, Filter_Name, model_pix_arcsec):
    """ Rebin an image for comparion to observations """
    if Filter_Name == 'HERSCHEL_PACS_BLUE':
        pix_arcsec = 1.6

    elif Filter_Name == 'HERSCHEL_PACS_RED':
        pix_arcsec = 3.2

    elif Filter_Name == 'JCMT_SCUBA2_450':
        pix_arcsec = 2.

    elif Filter_Name == 'JCMT_SCUBA2_850':
        pix_arcsec = 4.


    pix_area = pix_arcsec**2
    model_pix_area = model_pix_arcsec**2

    rebin_factor = pix_arcsec / model_pix_arcsec
    print(rebin_factor)
    rebinned_image = zoom(BeamConv_Image/model_pix_area, 1./rebin_factor)
    return rebinned_image
