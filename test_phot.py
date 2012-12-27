import PyGuide
import pyfits
import phot
import numpy
from process import flareCam
import time

file = './test.FIT'
img = pyfits.open(file)
data = img[0].data
img.close()

ccd = PyGuide.CCDInfo(0, flareCam['readNoise'], flareCam['ccdGain']) # from flarecam manual
refCoords, stats = PyGuide.findStars(data, None, None, ccd)

# look at 5 brightest
refCoords = numpy.asarray([refCoords[j].xyCtr for j in range(5)])
refCoords = refCoords - 0.5 # origin for imshow is effectively .5, .5 rather than 0,0


    
def experiment(inrad = 4, skyAnnulus = [7,9], 
                gridDense=31, splineOrder=0):
    photOut = []
    phot = phot.ApPhot(data, inrad = inrad, skyAnnulus = skyAnnulus, 
                gridDense=gridDense, splineOrder=splineOrder)
    for coord in refCoords:
        photOut.append(phot.fuckinDoIt(coord))
    return photOut
    
def loadTest(inrad = 4, skyAnnulus = [7,9], 
                gridDense=11, splineOrder=0):
    begTime = time.time()
    out = []
    for x in range(1000):
        if x % 10 == 0:
            print 'iter: ', x
        out.append(experiment(inrad = inrad, skyAnnulus = skyAnnulus, 
                gridDense=gridDense, splineOrder=splineOrder))
    print 'took %s for %i iterations' % (time.time()-begTime, x)
    return out
    
