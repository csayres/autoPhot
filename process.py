"""Do photometry

note: fix plotting shift
"""

import time
import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import PyGuide

from triangleHash import *
from preProcess import *
import phot

INRAD = 3 # pixels
SKYANN = [6,8] # pixels
RES = 21 # must be odd
SPLINE = 0

# flare cam header solutions, used in preprocessing
flareCam = {
'dateObs': 'DATE-OBS', # of form: 2012-05-20T11:29:23.338
'exptime': 'EXPTIME',
'filter': 'FILTER',
'imgExt': 'FIT', # image file extention
'readNoise': 11,
'ccdGain': 1.5 
}

# class StarNotFoundError(Exeception):
#     # implement later instead of runtimeErrors?
#     pass

class FieldPhot(object):
    """Photometry object for doing photometry on a field
    """
    def __init__(self, objName, dir = '../../data/', instDict = flareCam):
        """Inputs:
        dir - directory with data files to reduce.
        """
        self.objName = objName
        self.preProcess = PreProcess(objName, dir, instDict)
        self.targetCoords = None
        self.comparisonCoords = []
        self.ccd = PyGuide.CCDInfo(0, instDict['readNoise'], instDict['ccdGain']) # from flarecam manual
        self.refHash = None
        self.refData = None
        self.nBrightStars = 10 # only use this many stars for image alignments
        self.searchRad = None # pixels

    def circle(self, xyCtr, rad):
        """Takes a center and a radius, returns x and y for plotting
        """
        theta = numpy.linspace(0, 2*numpy.pi, 100)
        x = rad * numpy.cos(theta) + xyCtr[0]
        y = rad * numpy.sin(theta) + xyCtr[1]
        return x, y
        
    def showField(self, imgData, targ, comp, title = None, block = False, fig = None):
        """Show the ccd image
        
        imgData = 2D image array
        targ = None, or x,y position of a target to highlight
        ref = None, or list of x,y positions for all reference objects
        block = Plot blocks code
        """
        if fig !=None:
            plt.close(fig)
        fig = plt.figure()    
        axNew = fig.add_subplot(111)
        #norm = colors.Normalize()
        #axNew.imshow(imgData, vmin = 0, vmax = 1000, cmap=cm.Greys_r)
        axNew.imshow(numpy.log(imgData**6), cmap=cm.Greys_r)        
        # do circles
        targX, targY = self.circle(targ, self.searchRad)
        axNew.plot(targX, targY, 'r')
        for coords in comp:
            try:
                compX, compY = self.circle(coords, self.searchRad)
                axNew.plot(compX, compY, 'g')
            except:
                pass
            #axNew.plot((coords)[0], (coords)[1], 'og', ms = 20, alpha = 0.5)
        if title:
            plt.title(title)
        plt.show(block=block)
        time.sleep(0.5)
        return fig
        #plt.close(fig)

    def doPhot(self, showEvery = None):
        """Begin doing automatic photometry
        """
        if self.refHash == None:
            print 'No reference frame found. Must run setupField() to select stars first.'
            return   
        self.refHash = self.doHash(self.refdata)
        tstart = time.time()
        fig = plt.figure()
        for filter, objList in self.preProcess.filterDict.items():
            for num, obj in enumerate(objList):
                data = obj.data   
                # unnecessasry?             
                exptime = obj.exptime
                dateObs = obj.dateObs
                ######################
                newHash = self.doHash(data)
                targCent = self.centroid(self.targetCoords, data)
                try: 
                    if not targCent.isOK:
                        print 'target wasnt found'
                        raise RuntimeError('target not found')
                    offset = numpy.subtract(self.targetCoords, targCent.xyCtr) # if a small offset was present
                    for val in self.comparisonCoords:
                        cent = self.centroid(numpy.subtract(val, offset), data)
                        if not cent.isOK:
                            print 'comparison not found'
                            raise RuntimeError('comparison not found')# one of the comparison stars wasn't found, try solving for big offset
                except RuntimeError:                     
                    offset = self.refHash.hashItOut(newHash) # offset is computed every time
                    targCent = self.centroid(self.targetCoords - offset, data)
                    if not targCent.isOK:
                        print 'well Fuck, couldnt find target after offset'
                        fig = self.showField(data, self.targetCoords, self.comparisonCoords, block=True, fig = fig)
                        fig = self.showField(data, self.targetCoords-offset, self.comparisonCoords-offset, block=True, fig=fig)    
                        continue                    
                    # offset worked, update shit
                self.refHash = newHash
                self.targetCoords = targCent.xyCtr
                for ind, val in enumerate(self.comparisonCoords):
                    self.comparisonCoords[ind] = numpy.subtract(val, offset)
                compCent = []
                compShape = []
                compPhotObj = []
                apPhot = phot.ApPhot(data, inrad=INRAD, skyAnnulus = SKYANN, gridDense=RES, splineOrder = SPLINE)             
                for compCoord in self.comparisonCoords:
                    # do one reference star at a time
                    coord = self.centroid(compCoord, data)
                    if not coord.isOK:
                        print 'Comparison Warning: %s' % coord.msgStr
                        compPhotObj.append(None)
                    else:
                        compPhotObj.append(apPhot.fuckinDoIt(compCoord-0.5)) # PyGuide = 0.5 px shift
                    compCent.append(coord)
                    compShape.append(self.starShape(compCoord, data))
                if showEvery != None and num % showEvery == 0:      
                    fig = self.showField(data, self.targetCoords, self.comparisonCoords, title = 'Image Num %s' % num, fig = fig)
                # enter into dicionary
                self.preProcess.filterDict[filter][num].targCentroid = targCent
                self.preProcess.filterDict[filter][num].targShape = self.starShape(targCent.xyCtr, data)
                self.preProcess.filterDict[filter][num].compCentroid = compCent
                self.preProcess.filterDict[filter][num].compShape = compShape
                self.preProcess.filterDict[filter][num].compPhotObj = compPhotObj
                self.preProcess.filterDict[filter][num].targPhotObj = apPhot.fuckinDoIt(numpy.asarray(targCent.xyCtr)-0.5) #PyGuide = 0.5 px shift
                #self.showField(data, self.targetCoords, self.comparisonCoords)
        plt.close(fig)
        print 'took %s seconds' % (time.time() - tstart)
        
    def doHash(self, data):
        """Setup the triangle hash table for correcting image offsets
        Inputs:
        data: x,y image array
        """
        # Set up Hash Table for the 1st image (which all others will be compared to)
        refCoords, stats = PyGuide.findStars(data, None, None, self.ccd)
        #refCoords = numpy.asarray(refCoords) - 0.5 # origin is offset by half a pixel width
        # extract xy centers and just keep the nBrightStars
        num = min(len(refCoords), self.nBrightStars)
        refCoords = numpy.asarray([refCoords[j].xyCtr for j in range(num)])
        return TriangleHash(refCoords)
        
    def setupField(self, n=0, rad=None):
        """Select target and comparison stars
        Inputs:
        n = image number, defaults to zero
        """
        if rad != None:
            self.searchRad = rad
        # grab the first image from preprocessing
        # don't care about filter
        items = self.preProcess.filterDict.iteritems()
        filt, img = items.next()
        self.refdata = img[n].data
        self.refHash = self.doHash(self.refdata)
        # first select target star and comparisons
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #norm = colors.Normalize()
        ax.imshow(numpy.log(self.refdata**6), cmap=cm.Greys_r)
        #ax.imshow(self.refdata, vmin = 0, vmax = 1000, cmap=cm.Greys_r)
        ax.set_clim=(0.0,0.7)
        cid = fig.canvas.mpl_connect('button_press_event', self.onclick)
        print 'Click on target star'
        plt.show()
        if self.targetCoords == None:
            raise RuntimeError('No Target was selected!')
        if not self.comparisonCoords:
            print 'Warning: no comparision stars were selected'       

    def centroid(self, xyPos, data):
        """check to see if there is valid signal at xpos, ypos, on image data
        if so return the centroid.
        """
        return PyGuide.centroid(data, None, None, xyPos, rad = self.searchRad, #pixels
            ccdInfo = self.ccd)

    def starShape(self, xyPos, data):
        """fit a star shape at an xy position
        """
        return PyGuide.starShape(data, None, xyPos, self.searchRad)
        
    def onclick(self, event):
        """Click on plot callback function. Set the target star and comparision stars
        """
        cent = self.centroid([event.xdata, event.ydata], self.refdata)
        if cent.isOK:
            if self.targetCoords == None:
                # first click, it's the target
                self.targetCoords = numpy.asarray(cent.xyCtr)
                if self.searchRad == None: # set it automatically
                    self.searchRad = self.starShape(cent.xyCtr, self.refdata).fwhm * 3.
                #self.ax.plot(self.targCoords[0], self.targCoords[1], 'or', ms = 20, alpha = 0.5)
                print 'Got Target, now click on desired comarison stars.'
            else:
                # must be a comparison star click
                # check if we already got that target
                dist = numpy.linalg.norm(numpy.asarray(cent.xyCtr) - numpy.asarray(self.targetCoords))
                if dist < 4:
                    print 'Already got that star!'
                    return
                for comp in self.comparisonCoords:
                    dist = numpy.linalg.norm(numpy.asarray(cent.xyCtr) - numpy.asarray(comp))
                    if dist < 4:
                        print 'Already got that star!'
                        return
                self.comparisonCoords.append(numpy.asarray(cent.xyCtr))
                print 'Got comparison star'
        else:
            # centroid didn't work, try again?
            print 'Star not found...try another click?'
            
    def plotLightCurve(self, time, timeseries, title=None, xlabel=None, ylabel=None):
        """plot a light curve
        """
        plt.figure()
        timeseries = numpy.asarray(timeseries)
        plt.plot(time, timeseries, '.k')
        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)  
        if title:      
            plt.title(title)
        plt.ylim((.96, 1.06))
        plt.show()
        
    def plotTarg(self, filter = None):
        """Plot lightcurve for target
        """        
        if filter == None:
            # just use the first one...
            keys = self.preProcess.filterDict.keys()
            filter = keys[0]
        objs = self.preProcess.filterDict[filter]
        time = []
        timeseries = []
        for obj in objs:
            time.append(obj.dateObs)
            timeseries.append(obj.targPhot)
        time = numpy.asarray(time)
        timeseries = numpy.asarray(timeseries)
        self.plotLightCurve(
            time, timeseries, 
            title = 'Target Light Curve',
            xlabel = 'Date of Observation',
            ylabel = 'ADUs per second'
            )
    
    def plotRef(self, filter=None):
        """Plot a panel of all reference stars' lightcurves
        """
        if filter == None:
            # just use the first one...
            keys = self.preProcess.filterDict.keys()
            filter = keys[0]
        objs = self.preProcess.filterDict[filter]
        time = []
        timeseries = []
        for obj in objs:
            time.append(obj.dateObs)
            timeseries.append(obj.compPhot)
        time = numpy.asarray(time)
        timeseries = numpy.asarray(timeseries).T # n stars 0 dim, time 1 dim
        for num, ts in enumerate(timeseries):
            self.plotLightCurve(
                time, ts, 
                title = 'Target Light Curve %s' % str(num),
                xlabel = 'Date of Observation',
                ylabel = 'ADUs per second'
                )
                
    def plotDiff(self, filter=None):
        """Plot differential photometry for target
        """
        if filter == None:
            # just use the first one...
            keys = self.preProcess.filterDict.keys()
            filter = keys[0]
        objs = self.preProcess.filterDict[filter]
        time = []
        timeseries = []
        for obj in objs:
            time.append(obj.dateObs)
            timeseries.append(obj.diffPhot)
        time = numpy.asarray(time)
        timeseries = numpy.asarray(timeseries)
        m = numpy.median(timeseries)
        timeseries = timeseries / m # normalize so lightcurve sits at 1
        self.plotLightCurve(
            time, timeseries, 
            title = 'Target Light Curve',
            xlabel = 'Date of Observation',
            ylabel = '% Variation Flux'
            )

    def plotDiffPG(self, filter=None):
        """Plot differential photometry for target
        """
        if filter == None:
            # just use the first one...
            keys = self.preProcess.filterDict.keys()
            filter = keys[0]
        objs = self.preProcess.filterDict[filter]
        time = []
        timeseries = []
        for obj in objs:
            time.append(obj.dateObs)
            timeseries.append(obj.diffPhotPG)
        time = numpy.asarray(time)
        timeseries = numpy.asarray(timeseries)
        m = numpy.median(timeseries)
        timeseries = timeseries / m # normalize so lightcurve sits at 1
        self.plotLightCurve(
            time, timeseries, 
            title = 'Target Light Curve',
            xlabel = 'Date of Observation',
            ylabel = '% Variation Flux'
            )            
        

# for testing
# g = FieldPhot('test')
# g.setupField()
# g.doPhot(stack=True)
            
if __name__ == '__main__':
    g = FieldPhot('test')
    g.setupField()
    g.doPhot()
    #g.timeseries.plotGreen()
    #g.timeseries.plotDiff()