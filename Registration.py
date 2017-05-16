###################################################
## Functions and classes for image registration, ##
## processing and transformation                 ##
###################################################

import numpy as np
import skimage.io as skIO
import imreg_dft


def register_images(ref_image, query_image):
    '''
    Use the discrete Fourier transform to perform image registration.
    
    Scaling, rotation and translations are found in the frequency domain
    using the imreg_dft implementation.
    
    These scaling, rotation and translations are output to be applied to
    the matched hybridization images.
    
    Args
    ----
    ref_image: np.array
        Numpy array of pixel intensities for the reference image
        
    query_image: np.array
        Numpy array of pixel intensities for the query image
        
    Returns
    -------
    transform: dict
        dictionary containing keys `scale`:float, `angle`:float and `translate`:(x, y).
        Values are expressed in the time domain (original image co-ordinate space).
    '''
    
    # using imreg_dft for image registration
    # how many iterations are required? default is 3, is 1 enough?
    # need to check how robust/variable the transformation estimates are
    trans_dict = imreg_dft.imreg.similarity(ref_image, query_image, 1)
    scale = trans_dict['scale']
    rotation = trans_dict['angle']
    tvec = trans_dict['tvec']
    
    transform = {'scale': scale,
                 'angle': rotation,
                 'translate': tvec}
    
    return transform


def apply_transformation(hybe_image, scale=1.0, angle=0.0, x_trans=0.0, y_trans=0.0):
    '''
    Given some transformation (scale, rotation and translation), derived from
    the registration of fiduciary beads, apply this to the matching
    hybridizaiton image.
    
    If the `scale` paramater != 1, a warning is issued indicating the images
    should have been the same size during registration.
    
    If the `angle` parameter != 0 +/- eps (where -10 < eps < 10), then there
    may have been an issue during registration or image acquisition.
    
    Args
    ----
    hybe_image: np.array
        Numpy array of probe hybridisation signal for a given hybridisation
        
    scale: float
        Float, apply a scale factor to the image size (if not-1 then there mighy be issues with
        the registration)
        
    angle: float
        Rotation angle to apply.  MERFish manuscript indicates that non-zero rotations are due
        to failures in image registration.
        
    x_trans: float
        Translation to apply in the X co-ordinate space(rows)
        
    y_trans: float
        Translation to apply in the Y co-ordinate space (columns)
        
    Returns
    -------
    trans_image: np.array
        Numpy array of registered hybridisation image
    '''
    
    if not np.isclose(scale, 1.0, atol=0.001):
        logging.warning("The image scale parameter is not 1"
                        ", image sizes are different.  This may result "
                        "in unwanted image distortions")
    else:
        pass
    
    if not np.isclose(angle, 0.0, atol=0.5):
        logging.warning("The image rotation paramter, {}, is not 0"
                        ", there may have been an issue during fiduciary"
                        " bead registration or image acquisition".format(angle))
    else:
        pass
    
    trans_image = imreg_dft.imreg.transform_img(hybe_image, scale=scale,
                                                angle=angle, tvec=(x_trans, y_trans))
    
    return trans_image
