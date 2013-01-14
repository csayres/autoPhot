"""it's a test, ok?
"""

import pyfits
import numpy
import PyGuide

import smartPhot
from smartPhot.config import flareCam

ccdInfo = PyGuide.CCDInfo(0, flareCam.readNoise, flareCam.ccdGain)  
img = pyfits.open('/Users/csayres/data/arcsat/GJ1243/gFilter/oneNight/ccs20120520.00001139.LUMINANCE.FIT')
data = img[0].data
img.close()
coord1 = numpy.array([197.66827757,  209.33221882])
coord2 = numpy.array([ 343.39379969,  210.25893067])
cent1 = PyGuide.centroid(data, None, None, coord1, rad = 2, #pixels
            ccdInfo = ccdInfo)
cent2 = PyGuide.centroid(data, None, None, coord2, rad = 2, #pixels
            ccdInfo = ccdInfo)
diff1 = 

