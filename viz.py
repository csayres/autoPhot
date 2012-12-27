"""visualization tools
"""
import matploltlib.pyplot as plt

def apPhotExtraction(photObj):
    """Make a plot of the interpolated psf with overlays showing the radii of interest
    inputs:
    photObj: a phot.photObj object
    """
    #img = plt.imshow(photObj._img, interpolation='none', origin='upper')
    
    psfMask = numpy.zeros(photObj._img.shape)
    #imgOut = photObj._img[:]
    maskVal = numpy.median(photObj._counts) # so it scales to a nice color
    imgOut = numpy.zeros(photObj._img.shape) - maskVal
    skyMask = numpy.zeros(photObj._img.shape)
    
    for x,y in itertools.izip(photObj._countsInds[0], photObj._countsInds[1]):
        imgOut[x,y] = photObj._img[x,y]
    for x,y in itertools.izip(photObj._skyInds[0], photObj._skyInds[1]):
        imgOut[x,y] = photObj._img[x,y]
    img = plt.imshow(imgOut, interpolation='none', origin='upper')    
    #im2 = plt.imshow(psfMask, interpolation='none', origin='upper', alpha=.6)
    #im2 = plt.imshow(skyMask, interpolation='none', origin='upper', alpha=.6)   
    plt.plot(photObj._center[0], photObj._center[1], 'r.')
    plt.xlim(0, imgOut.shape[0]-1)
    plt.ylim(0, imgOut.shape[0]-1)     
    plt.colorbar()
    plt.show()
    
def skyHist(photObj):
    """Generate a histogram of the sky values
    inputs:
    photObj: a phot.photObj object
    """
    plt.hist(photObj._sky)
    plt.show()

def countsHist(fig, photObj):
    """Generate a histogram of the counts
    inputs:
    photObj: a phot.photObj object
    """
    plt.hist(photObj._counts)
    plt.show()
    

def plotDiff(fig, diffPhotObjs):
    """Plot differential photometry for target
    input:
    a list of DiffPhotObjs
    """
    time = []
    timeseries = []
    for obj in diffPhotObjs:
        time.append(obj.img.dateObs)
        timeseries.append(obj.diffCounts)
    time = numpy.asarray(time)
    timeseries = numpy.asarray(timeseries)
    m = numpy.median(timeseries)
    timeseries = timeseries / m # normalize so lightcurve sits around 1
    plotLightCurve(fig,
        time, timeseries, 
        title = 'Target Light Curve',
        xlabel = 'Date of Observation',
        ylabel = '% Variation Flux'
        )
        
def plotLightCurve(fig, time, timeseries, title=None, xlabel=None, ylabel=None):
    """plot a light curve
    """
    ax = fig.add_subplot(111)
    timeseries = numpy.asarray(timeseries)
    ax.plot(time, timeseries, '.k')
    if xlabel:
        ax.xlabel(xlabel)
    if ylabel:
        ax.ylabel(ylabel)  
    if title:      
        ax.title(title)
        
def showField(fig, diffPhotObj):
    """Show the ccd image
    
    inputs:
    diffPhotObj
    """
    ax = fig.add_subplot(111)
    #norm = colors.Normalize()
    ax.imshow(scaleImg(diffPhotObj.img.data), cmap=cm.Greys_r)        
    # do circles
    targX, targY = circle(diffPhotObj.targCentroid.xyCtr - 0.5, diffPhotObj.INRAD)
    ax.plot(targX, targY, 'r')
    for comp in diffPhotObj.compCentroids:
        try:
            compX, compY = self.circle(comp.xyCtr - 0.5, diffPhotObj.INRAD)
            axNew.plot(compX, compY, 'g')
        except:
            pass
        
def circle(xyCtr, rad):
    """Takes a center and a radius, returns x and y for plotting
    """
    theta = numpy.linspace(0, 2*numpy.pi, 100)
    x = rad * numpy.cos(theta) + xyCtr[0]
    y = rad * numpy.sin(theta) + xyCtr[1]
    return x, y

def scaleImg(imgData):
    """Figure out a good way to scale astro images?
    """
    return numpy.log(imgData**6)