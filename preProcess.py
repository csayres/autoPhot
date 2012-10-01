"""Look at files in the directory of interest.  Directory should contain
just one object.  Multiple filters are found and seperated.
"""

import pyfits
import collections
import itertools
import numpy
import datetime
import time
import glob
import matplotlib.pyplot
import matplotlib.dates


class ImageStackData(object):
    """Interface for holding data corresponding to a group of stacked images.
    """
    def __init__(self, imgDataList):
        """Inputs:
        imgDataList: a list of ImageData objects to be treated as a single (stacked)
        exposure.
        
        notes:
        date of observation is adopted from the first image in list
        filter is adopted from the first image in list
        exposure time is added for all exposures
        ccd data is simply stacked
        """
        exptime = 0
        self.fileList = []
        self.dateObs = imgDataList[0].dateObs
        self.filter =  imgDataList[0].filter
        for img in imgDataList:
            exptime += img.exptime
            self.fileList.append(img.filename)
        self.exptime = exptime
        self.targCentroid = None # to be entered upon photometry. PyGuide Centroid Object
        self.targShape = None # to be entered upon photometry. PyGuide StarShape Object
        self.compCentroid = [] # to be entered upon photometry. List of PyGuide Centroid Objects
        self.compShape = [] # to be entered upon photometry. List of PyGuide StarShape Objects
        
        
    @property
    def data(self):
        stackedData = None
        for file in self.fileList:
            img = pyfits.open(file)
            data = img[0].data
            img.close()
            if stackedData==None:
                stackedData = data
                continue
            stackedData += data
        return data

    @property        
    def compPhot(self):
        """Return background subtracted photometry for each comparison star
        output: list length(self.compCentorid), of ADU/second
        """ 
        photOut = []
        for cent, shape in itertools.izip(self.compCentroid, self.compShape):
            try:
                print 'ref pix: ', cent.pix
                photOut.append((cent.counts - (cent.imStats.med)*cent.pix)/self.exptime)
            except TypeError:
                photOut.append(numpy.nan)
        return numpy.asarray(photOut)
    
    @property        
    def targPhot(self):
        """Return background subtracted photometry for each comparison star
        output: a background subtracted ADU/second
        """ 
        cent = self.targCentroid
        print 'targ pix: ', cent.pix
        shape = self.targShape 
        return numpy.asarray((cent.counts - (cent.imStats.med)*cent.pix)/self.exptime)
        
    @property
    def diffPhot(self):
        """Return a differential flux/sec by dividing target flux by
        the sum of all reference star flux
        """
        refFluxAll = numpy.sum(self.compPhot)
        return self.targPhot / refFluxAll
                    
    
    
class ImageData(object):
    """an object for storing image attributes for easy lookup
    """
    def __init__(self, filename, dateObs, exptime, filter):
        """Inputs:
        filename = the filename
        dateObs = a python datetime object for the time of observation
        exptime = exposure time
        filter = filter
        """
        self.filename = filename
        self.dateObs = dateObs
        self.exptime = exptime
        self.filter = filter
        self.targCentroid = None # to be entered upon photometry. PyGuide Centroid Object
        self.targShape = None # to be entered upon photometry. PyGuide StarShape Object
        self.compCentroid = [] # to be entered upon photometry. List of PyGuide Centroid Objects
        self.compShape = [] # to be entered upon photometry. List of PyGuide StarShape Objects

    @property        
    def data(self):
        img = pyfits.open(self.filename)
        data = img[0].data
        img.close()
        return data

    @property        
    def compPhot(self):
        """Return background subtracted photometry for each comparison star
        output: list length(self.compCentorid), of ADU/second
        """ 
        photOut = []
        for cent, shape in itertools.izip(self.compCentroid, self.compShape):
            try:
                photOut.append((cent.counts - (cent.imStats.med)*cent.pix)/self.exptime)
            except TypeError:
                photOut.append(numpy.nan)
        return numpy.asarray(photOut)
    
    @property        
    def targPhot(self):
        """Return background subtracted photometry for each comparison star
        output: a background subtracted ADU/second
        """ 
        cent = self.targCentroid
        shape = self.targShape
        return numpy.asarray((cent.counts - (cent.imStats.med)*cent.pix)/self.exptime)
        
    @property
    def diffPhot(self):
        """Return a differential flux/sec by dividing target flux by
        the sum of all reference star flux
        """
        refFluxAll = numpy.sum(self.compPhot)
        return self.targPhot / refFluxAll

class PreProcess(object):
    """Go through all files and organize.
    """
    def __init__(self, objName, dir, instDict):
        """
        Inputs:
        objName = name of object
        dir = directory where object fits files are stored
        instDict = an insturment dictory containing 
            necessary info written by a given isntrument
        """
        self.objName = objName
        allObjs = []
        filters = []
        # grab every image file in the directory
        allFitsFiles = glob.glob(dir + '*.' + instDict['imgExt'])
        for file in allFitsFiles:
            # extract and save useful header data from each image
            imgFile = pyfits.open(file)
            filter = imgFile[0].header[instDict['filter']].strip() # get filter used, strip whitespace
            dateObs = self.convertDateStr(imgFile[0].header[instDict['dateObs']])
            exptime = imgFile[0].header[instDict['exptime']]
            if filter not in filters:
                filters.append(filter)
            allObjs.append(
                ImageData(file, dateObs, exptime, filter)               
            )
            imgFile.close()
        # organize this shit in a useful way: by filter and observation date
        self.filterDict = self.sort2Dict(allObjs, filters)

    def sort2Dict(self, objs, filters):
        """Sort objects into a dicionary keyed by each filter, 
        Input:
        objs - a list of ImageData objects
        filters - filters to key the dictionary with
        """
        # sort by datetime filter and by datetime of observaiton
        order = numpy.argsort([obj.dateObs for obj in objs])
        allObjs = [objs[a] for a in order] # reorder by time of observation
        filterDict = {}
        for filter in filters:
            filterDict[filter] = []
        for obj in allObjs:
            filterDict[obj.filter].append(obj)
        return filterDict

    def newDictStacked(self, oldDict, stack = 'auto'):
        """Take a current filter keyed object dictionary, 
        and apply stacking of your choice to it.
        
        oldDict: a dictionary keyed by filter names, of lists of ImageData
            objects in chronological order
        stackType: executable, one of (self.autoGroupImgs or self.numGroupImgs)
        stack: 'auto' or an integer.  If auto, and automatic stacking solution
            will be found. Or provide the number of frames to stack
        """
        newDict = {}
        for filt, imgs in oldDict.iteritems():
            if stack == 'auto':
                newDict[filt] = self.autoGroupImgs(imgs)
            else:
                newDict[filt] = self.numGroupImgs(stack)
        return newDict
            
    def getImgData(self, stack = None):
        """Return a dictionary of:
        'filter': List of ImageStackData (stacked) or ImageData objects
        """
        if stack != None:
            try:
                stack = int(stack)
                return self.newDictStacked(self.filterDict, stack)
            except ValueError:
                if stack != 'auto':
                    raise RuntimeError('stack value must be an integer, "auto", or None')
                return self.newDictStacked(self.filterDict, stack)
        else:
            # no stacking asked for, return every object
            return self.filterDict     
     
    def autoGroupImgs(self, imgDataList):
        """Will return a list of subarrays of imgDataList elements.  Subsequent
        exposures will be determined based on exposure time gaps > some number.
        Should be useful when some sort of cycling through filters was used.
        
        Input:
        imgDataList: a list of ImageData objects sorted by exposure time
        
        Output:
        a list of ImageStackData objects
        """
        multiplier = 2.5 # throw out exposures that are seperated by more than this factor of exptime        
        exptime = numpy.median(numpy.asarray([img.exptime for img in imgDataList]))
        obstimes = numpy.asarray([img.dateObs for img in imgDataList])
        diffs = numpy.asarray([x.total_seconds() for x in numpy.diff(obstimes)])
        cutoff = numpy.median(diffs) * multiplier  
        # nonzero returns a tuple..hence the [0] index to get the array      
        splits = numpy.nonzero(diffs > cutoff)[0] + 1
        groupList = numpy.split(numpy.asarray(imgDataList), splits)
        return [ImageStackData(x) for x in groupList]
        
    def numGroupImgs(self, imgDataList, num):
        """return imgDatalist broken into sub-arrays of lengh num
        if imgDataList isn't divisible by num, the remainder of images will
        be tossed.
        
        Input:
        imgDataList: a list of ImageData objects sorted by exposure time
        num: a number of images to return in each subarray
        
        Output
        a list of ImageStackData objects
        """
        rem = len(imgDataList) % num
        if  rem !=0:
            #throw away the images in the remainder
            n = len(imgDataList[-1*rem:])
            del imgDataList[-1*rem:]
            print 'warning: throwing away %s leftover images due to stacking remainder' % n
        groupList = numpy.split(numpy.asarray(imgDataList), num)
        return [ImageStackData(x) for x in groupList]
        
    def convertDateStr(self, dateStr):
        """Parse and convert a date string from the flare-cam header 
        to a python datetime object
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

if __name__ == '__main__':        
    test = PreProcess('test')
            