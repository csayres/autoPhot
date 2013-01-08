"""visualization tools
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy
import itertools

# axNew = fig.add_subplot(111)

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
    

def plotDiff(ax, diffPhotObjs):
    """Plot differential photometry for target
    input:
    ax: from matplotlib.pyplot.figure.add_subplot(...)
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
    plotLightCurve(ax,
        time, timeseries, 
        title = 'Target Light Curve',
        xlabel = 'Date of Observation',
        ylabel = '% Variation Flux'
        )
        
def plotLightCurve(ax, time, timeseries, title=None, xlabel=None, ylabel=None):
    """plot a light curve
    """
    ax.plot(time, timeseries, '.k')
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)  
    if title:      
        plt.title(title)
        
def showFieldSolution(ax, fieldSolution):
    """display field solution data, must input a figure and FieldSolution object
    """
    circRad = 5
    #norm = colors.Normalize()
    ax.imshow(scaleImg(fieldSolution.img.data), cmap=cm.Greys_r)
    if fieldSolution.targetCoords != None:
        targX, targY = circle(numpy.asarray(fieldSolution.targetCoords) - 0.5, circRad)
        ax.plot(targX, targY, 'r')
    if fieldSolution.compCoords != None:
        for comp in fieldSolution.compCoords:
            compX, compY = circle(numpy.asarray(comp) - 0.5, circRad)
            ax.plot(compX, compY, 'g')
    plt.xlim((0,fieldSolution.img.data.shape[0]))   
    plt.ylim((0,fieldSolution.img.data.shape[1]))        
    plt.draw()
        
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