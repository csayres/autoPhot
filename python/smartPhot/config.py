"""Configuration elements for photometry, etc
"""
import datetime
import numpy
import img
import PyGuide

class PhotConfig(object):
    """Simple object that holds configuration for photometry
    
    includes defaults for photometry
    """
    def __init__(self,           
            inrad = 2, # pixels
            skyann = [6, 8], # pixels
            res=21, # must be odd
            spline=0, # should probably remain 0.  Higher orders seem to
            # introduce some structure (aliasing?) into the signal
            # more experimentation is probably merited.
        ):
        self.inrad = inrad
        self.skyann = skyann
        self.res = res
        self.spline = spline

class CameraConst(object):
    """A dictionary for storing camera attributes
    """
    def __init__(self, dateObs, exptime, filter, imgExt, readNoise, ccdGain):
        """Inputs (strings corresponding to image header):
        """
        self.dateObs = dateObs
        self.exptime = exptime
        self.filter = filter
        self.imgExt = imgExt
        self.readNoise = readNoise
        self.ccdGain = ccdGain
    
    def parseDate(self, dateStr):
        """defined by subclasses, must return a datetime object, constructed by 
        parsing the dateObs string from the image header specific to the 
        camera.
        """
        raise NotImplementedError('subclasses must override')

def parseDate(dateStr):
    """Parse and convert a date string from the flare-cam header 
    to a python datetime object
    of form: 2012-05-20T11:29:23.338
    """
    date, time = dateStr.split('T')
    year, month, day = date.split('-')
    # round to the nearest second
    hour, minute, second = time.split(':')
    year = int(year)
    month = int(month)
    day = int(day)
    hour = int(hour)
    minute = int(minute)
    second = int(numpy.round(float(second)))
    # seconds must be between 0 and 59 for datetime...
    if second == 60:
        second = 0
        minute += 1       
    return datetime.datetime(year, month, day, hour, minute, second) 

flareCam = CameraConst(
        dateObs = 'DATE-OBS', 
        exptime = 'EXPTIME', 
        filter = 'FILTER', 
        imgExt = 'FIT', 
        readNoise = 11, 
        ccdGain = 1.5,
        )
flareCam.parseDate = parseDate

class Config(object):
    """object containing all the necessary parameters for photometry, etc
    """
    def __init__(self, camera, phot, objFileList, biasFileList = None, flatFileList = None):
        """Inputs:
        camera: a CameraConst object
        phot: a PhotConfig object
        objFileList = list of files to process
        biasFileList = list of bias image files or None
        flatFileList = list of flat image files or None
        """
        self.camera = camera
        self.phot = phot  
        calibrator = None
        if biasFileList and flatFileList:
            biasList = img.imgLister(biasFileList, type = 'bias', 
                                        cameraConst = camera, calibrator = None)
            flatList = img.imgLister(flatFileList, type = 'flat', 
                                cameraConst = camera, calibrator = None)
            calibrator = img.Calibrator(biasList, flatList) # removes bias from flats
        else:
            print 'no calibration frames reveived!, proceeding'
        self.calibrator = calibrator
        self.objList = img.imgLister(objFileList, type = 'light', 
                            cameraConst = camera, calibrator = calibrator)
        self.ccdInfo = PyGuide.CCDInfo(0, camera.readNoise, camera.ccdGain)    
          
#FlareCamConfig = 