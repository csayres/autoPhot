"""Look at files in the directory of interest.  Directory should contain
just one object.  Multiple filters are found and seperated.

note: pyguide set to 5 px search radius
note: add way to check if reference is outside image to avoid unnecessary hashing.
note: add warning if comparison == target
note: to deal with shift during exposure, grab oldest field solution
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
#from config import flareCamConfig
import viz

from triangleHash import *

NBRIGHTSTARS = 20

class Driver(object):
    """Drives the reduction process, keeps track of what's going on and what's next
    """
    def __init__(self, config): #objFileList, config = flareCamConfig):#, cameraConst = flareCam, biasFileList = None, flatFileList = None):
        """biasList: list of bias images to use, a path
        flatList: list of flats to use, a path
        objPathList: list of obj images to use, should all be of same field, path
        """
        self.config = config
        self.fieldSolution = FieldSolution() # initialize empty
        self.cruncher = None # set by crunch method
        self.df = None # set by crunch method
        self.cruncher = Cruncher(self.fieldSolution, self.config.ccdInfo)
        self.cruncher.centroid = self.centroid

    def crunch(self):
        self.df = self.cruncher.crunchLoop(self.config.objList, self.config)
        
    def chooseTarget(self, imgNum = 0):
        """Select the target star
        """
        img = self.config.objList[imgNum]
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
            
        self.fieldSolution.img = img
        self.fieldSolution.hash = TriangleHash(self.cruncher.findSources(img.data))
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
            
        self.fieldSolution.img = self.config.objList[imgNum]
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
        return PyGuide.centroid(data, None, None, xyPos, rad = 3, #pixels
            ccdInfo = self.config.ccdInfo)
        
class Cruncher(object):
    """This does all the work, applying calibrations, finding stars, doing photometry...
    """
    def __init__(self, fieldSolution, ccd):
        """Inputs:
        fieldSolution: a FieldSolution object, which contains a triangleHash and locations for 
            target and comparison stars. Will get updated upon iterations
        imgList: a list of img.Light objects for photometry extraction
        ccd: a PyGuide ccd constant object
        """ 
        self.fieldSolution = fieldSolution
        self.ccd = ccd      
                
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
        
    def crunchPhot(self, img, offset, config):
#     def crunchPhot(self, img, config):
        """Do photometery for a single image
        returns a diffPhot Obj and an updated fieldSolution
        returns None if there was some problem
        """
        imgData = img.data
        
        ##########new#############
#         newSources = self.findSources(imgData)
#         try:            
#             offset = numpy.asarray(self.fieldSolution.hash.hashItOut(newSources))
#         except Exception as e:
#             print 'offset failed with exception %s'%e
#             return None
#         print 'new offset: ', offset
#         fig = plt.gcf()
#         #fig = plt.figure()
#         ax = fig.add_subplot(111)
#         viz.showField(ax, imgData, self.fieldSolution.targetCoords - offset, [comp - offset for comp in self.fieldSolution.compCoords])
#         plt.show(block=False)
        #############################        
        # check to see if all centroids are found in new image
        targCent = self.centroid(numpy.asarray(self.fieldSolution.targetCoords) - offset, imgData)
        compCent = [self.centroid(numpy.asarray(comp) - offset, imgData) for comp in self.fieldSolution.compCoords]  
        allCents =  compCent[:]
        allCents.append(targCent)  
        if False in [cent.isOK for cent in allCents]:
            ####new####
#             print 'sources missing'
#             return None  
            ###########
            # one or more centroids failed, try to hash for an offset
            newSources = self.findSources(imgData)
            tempFieldSolution = FieldSolution(
                targCent.xyCtr, 
                [comp.xyCtr for comp in compCent],
                img,
                )
            # overwrite the old offset
            try:            
                offset = numpy.asarray(self.fieldSolution.hash.hashItOut(newSources))
            except Exception as e:
                print 'offset failed with exception %s'%e
                return None
            print 'new offset: ', offset   
            targCent = self.centroid(self.fieldSolution.targetCoords - offset, imgData)
            compCent = [self.centroid(comp - offset, imgData) for comp in self.fieldSolution.compCoords]
            fig = plt.gcf()
            #fig = plt.figure()
            ax = fig.add_subplot(111)
            viz.showField(ax, imgData, self.fieldSolution.targetCoords - offset, [comp - offset for comp in self.fieldSolution.compCoords])
            plt.show(block=False)
            if not targCent.isOK:
                print "couldn't find target in image: %s" % img.path
                return None
            if False in [comp.isOK for comp in compCent]:
                print "warning: one or more comparison stars not found in image: %s" % img.path
                return None
            # target was found, keep the offset and run with it                
        try:
            df = phot.DiffPhotObj(
                img = img, targCentroid = targCent, 
                compCentroids = compCent, inrad = config.phot.inrad, skyAnnulus = config.phot.skyann, 
                resolution = config.phot.res, spline = config.phot.spline
                )
            return df, offset        
        except Exception as e:
            print 'Could not compute photometry, exception: %s' %e    
            return None 

    def crunchLoop(self, imgList, config):
        """begin the crunching, build a list of DiffPhot objects
        
        should be parrallelizable
        """
        # not sure if this is necessary, rationale is for parallel loops 
        # each with their own independent fieldSolution 
        #fieldSolution = copy.copy(self.fieldSolution) 
        dfList = []
        offset = numpy.array([0,0]) # begin with no offset
        imNum = 1
        for img in imgList:
                print 'image Number: ', imNum
                imNum += 1
                out = self.crunchPhot(img, offset, config)
                #out = self.crunchPhot(img, config)
                if not out:
                    # photometry returned None, skip that image
                    # try again using original field solution
                    print 'image extraction failed %s, skipping' % img.path
                    continue
                # update offset
                # append the diffPhotObj to the list
                df, offset = out
                dfList.append(df)
        return dfList
 
    def findSources(self, imgData):
        """get a list of automatically detected sources
        """ 
        # Set up Hash Table for the 1st image (which all others will be compared to)
        refCoords, stats = PyGuide.findStars(imgData, None, None, self.ccd)
        #refCoords = numpy.asarray(refCoords) - 0.5 # origin is offset by half a pixel width
        # extract xy centers and just keep the NBRIGHTSTARS
        num = min(len(refCoords), NBRIGHTSTARS)
        refCoords = numpy.asarray([refCoords[j].xyCtr for j in range(num)])
        return refCoords
    
    def computeOffset(self, imgData):
        """needed?
        """ 
        coordList = self.findSources(imgData)
        self.fieldSolution.hash.hashItOut(coordList)
        
class FieldSolution(object):
    """Contains information about where to find target and reference stars, and a hash
    table for the image in case of shifts
    """
    def __init__(self, targetCoords = None, compCoords = None, img = None, hash = None):
        """inputs:
        triangleHash: a TriangleHash object, used to determine image shifts
        targetCoords: [x,y] coordset
        compCoords: 2xN array of comparison coordinates
        img: the img.imgBase image used
        """
        self.targetCoords = targetCoords
        self.compCoords = compCoords
        self.img = img
        self.hash = hash 
            