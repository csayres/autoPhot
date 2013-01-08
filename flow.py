"""Look at files in the directory of interest.  Directory should contain
just one object.  Multiple filters are found and seperated.

note: pyguide set to 5 px search radius
note: add way to check if reference is outside image to avoid unnecessary hashing.
note: add warning if comparison == target
"""

import pyfits
import numpy
import datetime
import glob
import matplotlib.pyplot as plt
plt.ion()
import matplotlib.dates
import PyGuide
import copy
import phot
import img
from camera import flareCam
import viz

from triangleHash import *

INRAD = 3 # pixels
SKYANN = [6,8] # pixels
RES = 21 # must be odd
SPLINE = 0 # should probably remain 0.  Higher orders seem to
# introduce some structure (aliasing?) into the signal
# more experimentation is probably merited.
NBRIGHTSTARS = 10

class Driver(object):
    """Drives the reduction process, keeps track of what's going on and what's next
    """
    def __init__(self, objFileList, cameraConst = flareCam, biasFileList = None, flatFileList = None):
        """biasList: list of bias images to use, a path
        flatList: list of flats to use, a path
        objPathList: list of obj images to use, should all be of same field, path
        """
        calibrator = None
        if biasFileList and flatFileList:
            biasList = img.imgLister(biasFileList, type = 'bias', 
                                        cameraConst = flareCam, calibrator = None)
            flatList = img.imgLister(flatFileList, type = 'flat', 
                                cameraConst = flareCam, calibrator = None)
            calibrator = img.Calibrator(biasList, flatList) # removes bias from flats
        else:
            print 'no calibration frames reveived!, proceeding'
        self.calibrator = calibrator
        self.objList = img.imgLister(objFileList, type = 'light', 
                            cameraConst = flareCam, calibrator = calibrator)
        self.ccdInfo = PyGuide.CCDInfo(0, cameraConst.readNoise, cameraConst.ccdGain) # from flarecam manual
        self.fieldSolution = FieldSolution() # initialize empty
        self.cruncher = None # set by crunch method
        self.df = None # set by crunch method
        self.cruncher = Cruncher(self.fieldSolution, self.objList, self.ccdInfo)
        self.cruncher.centroid = self.centroid

    def crunch(self):
        self.df = self.cruncher.crunchLoop(self.objList)
        
    def chooseTarget(self, imgNum = 0):
        """Select the target star
        """
        def setTarget(event):
            """matplotlib clicky callback
            """
            cent = self.centroid(numpy.asarray([event.xdata, event.ydata]) + 0.5, self.fieldSolution.img.data)
            if cent.isOK:
                self.fieldSolution.targetCoords = numpy.asarray(cent.xyCtr)
                print 'got target'
                viz.showFieldSolution(event.inaxes, self.fieldSolution)
            else:
                print 'target not found, try again?'
            
        self.fieldSolution.img = self.objList[imgNum]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        viz.showFieldSolution(ax, self.fieldSolution)
        cid = fig.canvas.mpl_connect('button_press_event', setTarget)
        print 'Click on target star'
            
    def chooseComparisons(self, imgNum = 0):
        """Select the target star
        """
        compCoords = []
        def setComp(event):
            """matplotlib clicky callback
            """
            cent = self.centroid(numpy.asarray([event.xdata, event.ydata]) + 0.5, self.fieldSolution.img.data)
            if cent.isOK:
                compCoords.append(numpy.asarray(cent.xyCtr))
                self.fieldSolution.compCoords = compCoords # overwrite each time
                print 'got comparison star'
                viz.showFieldSolution(event.inaxes, self.fieldSolution)
            else:
                print 'comparison not found, try again?'
            #self.fieldSolution.compCoords = compCoords if compCoords else None
            
        self.fieldSolution.img = self.objList[imgNum]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        viz.showFieldSolution(ax, self.fieldSolution)
        cid = fig.canvas.mpl_connect('button_press_event', setComp)
        print 'Click on comparison star'
        #plt.show()
        #self.fieldSolution.compCoords = compCoords if compCoords else None

    def centroid(self, xyPos, data):
        """check to see if there is valid signal at xpos, ypos, on image data
        if so return the centroid.
        """
        # default to 5 px search rad?
        return PyGuide.centroid(data, None, None, xyPos, rad = 8, #pixels
            ccdInfo = self.ccdInfo)
        
class Cruncher(object):
    """This does all the work, applying calibrations, finding stars, doing photometry...
    """
    def __init__(self, fieldSolution, imgList, ccd):
        """Inputs:
        fieldSolution: a FieldSolution object, which contains a triangleHash and locations for 
            target and comparison stars. Will get updated upon iterations
        imgList: a list of img.Light objects for photometry extraction
        ccd: a PyGuide ccd constant object
        """ 
        self.fieldSolution = fieldSolution
        self.imgList = imgList
        self.ccd = ccd
        self.setPhotOptions()

    def setPhotOptions(self,            
            inrad = INRAD, # pixels
            skyann = SKYANN, # pixels
            res=RES, # must be odd
            spline=SPLINE, # should probably remain 0.  Higher orders seem to
            # introduce some structure (aliasing?) into the signal
            # more experimentation is probably merited.
            nbrightstars=NBRIGHTSTARS,
        ):
        self.inrad = inrad
        self.skyann = skyann
        self.res = res
        self.spline = spline
        self.nbrightstars = nbrightstars         
                
    def centroid(self, xyPOs, data):
        """Added by the Driver, 
        """
        raise NotImplementedError

    def inImgBound(self, coord, imgData):
        """Check that coord (xyPos) lies within the image bounds. 
        Returns True or False
        take PyGuide shift into account
        """
        coord = coord - 0.5
        return True if numpy.zeros(2) <= coord <= coord.shape else False
        
    def crunchPhot(self, img, fieldSolution):
        """Do photometery for a single image
        returns a diffPhot Obj and an updated fieldSolution
        returns None if there was some problem
        """
        imgData = img.data
        newHash = self.doHash(imgData) # might be able to get rid of this?
        # check to see if all centroids are found in new image
        targCent = self.centroid(fieldSolution.targetCoords, imgData)
        compCent = [self.centroid(comp, imgData) for comp in fieldSolution.compCoords]  
        allCents =  compCent[:]
        allCents.append(targCent)  
        if False in [cent.isOK for cent in allCents]:
            fig = plt.gcf()
            #fig = plt.figure()
            ax = fig.add_subplot(111)
            viz.showFieldSolution(ax, fieldSolution)
            plt.show(block=False)
            # one or more centroids failed, try to hash for an offset
            offset = fieldSolution.hash.hashItOut(newHash)    
            print 'new offset: ', offset   
            targCent = self.centroid(fieldSolution.targetCoords - offset, imgData)
            compCent = [self.centroid(comp - offset, imgData) for comp in fieldSolution.compCoords]
            if not targCent.isOK:
                print "couldn't find target in image: %s, ignoring" % img.path
                return None
            if False in [comp.isOK for comp in compCent]:
                print "warning: one or more comparison stars not found in image: %s" % img.path
            # target was found, keep the offset and run with it                

        df = phot.DiffPhotObj(img = img, targCentroid = targCent, 
            compCentroids = compCent, inrad = self.inrad, skyAnnulus = self.skyann, 
            resolution = self.res, spline = self.spline
        )
        newFieldSolution = FieldSolution(newHash, 
            targCent.xyCtr, 
            [comp.xyCtr for comp in compCent],
            img,
        )
        return df, newFieldSolution               

    def crunchLoop(self, imgList):
        """begin the crunching, build a list of DiffPhot objects
        
        should be parrallelizable
        """
        # not sure if this is necessary, rationale is for parallel loops 
        # each with their own independent fieldSolution 
        fieldSolution = copy.copy(self.fieldSolution) 
        dfList = []
        for img in imgList:
            out = self.crunchPhot(img, fieldSolution)
            if not out:
                # photometry returned None, skip that image
                continue
            else:
                # update FieldSolution
                # append the diffPhotObj to the list
                df, fieldSolution = out
                dfList.append(df)
        return dfList
                

    def doHash(self, data):
        """create a triangle hash from data
        inputs:
        data: 2d image array
        """
        # Set up Hash Table for the 1st image (which all others will be compared to)
        refCoords, stats = PyGuide.findStars(data, None, None, self.ccd)
        #refCoords = numpy.asarray(refCoords) - 0.5 # origin is offset by half a pixel width
        # extract xy centers and just keep the NBRIGHTSTARS
        num = min(len(refCoords), self.nbrightstars)
        refCoords = numpy.asarray([refCoords[j].xyCtr for j in range(num)])
        return TriangleHash(refCoords)    


class FieldSolution(object):
    """Contains information about where to find target and reference stars, and a hash
    table for the image in case of shifts
    """
    def __init__(self, triangleHash = None, targetCoords = None, compCoords = None, img = None):
        """inputs:
        triangleHash: a TriangleHash object, used to determine image shifts
        targetCoords: [x,y] coordset
        compCoords: 2xN array of comparison coordinates
        img: the img.imgBase image used
        """
        self.hash = triangleHash
        self.targetCoords = targetCoords
        self.compCoords = compCoords
        self.img = img

class PreProcess(object):
    """Go through all files and organize.
    """
    def __init__(self, objName, dir, cameraConst):
        """
        Inputs:
        objName = name of object
        dir = directory where object fits files are stored
        cameraConst = an insturment dictory containing 
            necessary info written by a given isntrument
        """
        self.objName = objName
        allObjs = []
        filters = []
        # grab every image file in the directory
        allFitsFiles = glob.glob(dir + '*.' + cameraConst.imgExt)
        for file in allFitsFiles:
            # extract and save useful header data from each image
            imgFile = pyfits.open(file)
            filter = imgFile[0].header[cameraConst.filter.strip()] # get filter used, strip whitespace
            dateObs = self.convertDateStr(imgFile[0].header[cameraConst.dateObs])
            exptime = imgFile[0].header[cameraConst.exptime]
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
            