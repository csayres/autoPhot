"""camera specifics
"""
import datetime
import numpy

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
    

class FlareCam(CameraConst):
    def __init__(self, 
        dateObs = 'DATE-OBS', 
        exptime = 'EXPTIME', 
        filter = 'FILTER', 
        imgExt = 'FIT', 
        readNoise = 11, 
        ccdGain = 1.5,
    ):
        CameraConst.__init__(self, dateObs, exptime, filter, imgExt, readNoise, ccdGain)

    def parseDate(self, dateStr):
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

flareCam = FlareCam()