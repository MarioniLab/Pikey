###################################################
## Dirty wrapper for MERFish example data files  ##
## so that they can be smushed into scikit-image ##
###################################################

import skimage
import skimage.io as skIO
import os
import sys

# find files
# use skimage for file IO and processing
# assume raw images in standard format, i.e. tiff, png, etc, not a custom/proprietary format
# images are 256x256 as standard?
# What class/object should be used for handling processing post-image acquisition?
# Operate entirely on a single object, or separate tasks to aid parallelisation?
# Can Python objects be serialized to HDF5?
# how much image format parsing needs to be done?
