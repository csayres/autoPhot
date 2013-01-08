"""Tools for doing aperture photometry

note: what do do if center is not inside image?

note: PyGuide introduces a 0.5 pixel shift on centroiding (because pixels are
defined with centers at n*(1-0.5) rather than n

photometry set to nan if pyguide centroid is not ok--don't think it will ever
be not ok, but whatever.
"""

import numpy
import scipy.ndimage as nd
import matplotlib.pyplot as plt
import itertools

class PhotObj(object):
    """A class for holding the photometry information
    """
    def __init__(self, countsInfo, skyInfo, img, center, radii, gridDense):
        """Inputs:
        countsInfo: [counts, inds] vector containing all counts for partial pixels, inds 2d array with xy indices, shold be background subtracted
        skyInfo: [counts, inds] vector containing all sky values for partial pixels, inds 2d array with xy indices
        img: the interpolated image region used for this specific psf
        center: the determined center of the psf, scaled to match img
        radii: [aperture, innerSky, outerSky], in original pixel scaling
        gridDense: the gridDense factor used to create partial pixels
        """
        self._counts = countsInfo[0]
        self._countsInds = countsInfo[1]
        self._sky = skyInfo[0]
        self._skyInds = skyInfo[1]
        self._center = center
        self._img = img
        self._radii = radii
        self._gridDense = gridDense

    @property
    def counts(self):
        """Return the number of counts for the psf, summing partial pixel values.
        """
        return numpy.sum(self._counts)

    @property
    def sky(self):
        """Return the sky value, this is rescaled (backwards) to match level corresponding
        to an original pixel (not a partial pixel).
        """
        return self.skyPartial * (self._gridDense**2)
    
    @property
    def skyPartial(self):
        """median sky value of a partial pixel
        """         
        return numpy.median(self._sky)
    

class ApPhot(object):
    """Object for doing aperture photometry
    """
    def __init__(self, data, inrad, skyAnnulus, 
                gridDense=11, splineOrder=0):
        """
        data: 2d image array
        center: xy center of object (psf)
        inrad: inner radius aperture
        skyAnnulus: (inner radius, outer radius)
        gridDense: increase grid density by this factor, determines number of partial
            pixels to be made
        splineOrder: order of the spline to fit to the upsampled image section.
                    0 = none.    
        """
        self.data = data
        if gridDense % 2 == 0:
            gridDense += 1 # I told you it must be odd.
        self.inrad = inrad 
        self.skyAnnulus = skyAnnulus 
        self.gridDense = gridDense
        self.splineOrder = splineOrder
    
    def fuckinDoIt(self, center):
        """Do aperture photmetry around a center point, return a photObj
        """
        center = numpy.asarray(center)
        # keep only region of image centered around object of square length about
        # the outer annulus
        region = self._regionExtract(center)
        # create a denser version for partial pixel emulation
        denseData = self._denseify(region)
        # center in new coordinates
        denseCenter = self._getNewCent(center, denseData)
        # get sky photometry
        # Annulus was given in pixels corresponding to the original image data
        # since then we have extracted a region and densified it, so the Annulus
        # values must be scaled accordingly
        skyInfo = self.radialExtract(denseData, denseCenter, 
                                            self.skyAnnulus[1] * self.gridDense,
                                            self.skyAnnulus[0] * self.gridDense)
        # determine median sky value and subtract
        skyLevel = numpy.median(skyInfo[0])
        denseNoSky = denseData - skyLevel
        countsInfo = self.radialExtract(denseNoSky, denseCenter, 
                                        self.inrad * self.gridDense)
        return PhotObj(countsInfo = countsInfo, 
                        skyInfo = skyInfo, 
                        img = denseNoSky, 
                        center = denseCenter, 
                        radii = [self.inrad, self.skyAnnulus[0], self.skyAnnulus[1]], 
                        gridDense = self.gridDense)
        
        

    def _getNewCent(self, oldCenter, interpData):
        """Determine the new center position of the psf after regionExtract/densify
        Inputs:
        center: the originial xyPos of the psf corresponding to self.data
        interpData: the new region-extracted and densified data
        
        Returns:
        A new center point corresponding to interpData
        """
        # offset of center from it's nearest pixel (which was used as the center
        # pixel for _regionExtract, which then must be scaled by the new grid density
        scaledOffset = (numpy.round(oldCenter) - oldCenter)*self.gridDense
        # new centerpoint is the coords of this offset from the center of interpData
        return numpy.floor(interpData.shape[0]/2.) - scaledOffset       
        
        
    def _regionExtract(self, center):
        """Extract square region around star of interest. Box-size 
        based on skyAnnulus.
        
        returns
        region: the extracted square region, centered on the nearest pixel
            to center. Sidelength is odd and equal to self.skyAnnulus[1] or
            self.skyAnnulus[1]+1 (whichever is odd)
        """
        boxSize = self.skyAnnulus[1]
        # box must have odd size length, so there is always 1 center pixel
        if self.skyAnnulus[1] % 2 == 0:
            boxSize = boxSize + 1
        centerPix = numpy.round(center) # pixel closest to psf center
        # get region of image centered on centerPix, of size (boxSize, boxSize)
        # note x,y reversal, it works.
        region = self.data[centerPix[1]-boxSize:centerPix[1]+boxSize+1, centerPix[0]-boxSize:centerPix[0]+boxSize+1]
        return region
        
    def _denseify(self, region):
        """slice the image into a denser grid. Should be done after _regionExtract,
        so you're working with less pixels, and only the region of interest.
        
        if self.spline order > 0, smoothing will occur during this process
        
        inputs:
        region: the image data, should be a region rather than a whole image
        
        returns: 
        denseData: the upsampled (and possibly interpolated) data 
        """ 
        denseData = nd.zoom(region, self.gridDense, order=self.splineOrder, prefilter=False)
        # grid is higher density so must renormalize so each pixel counts for less
        denseData = denseData/(self.gridDense**2) 
        return denseData
                    
    def radialExtract(self, data, center, outer, inner=0):
        """Return values from data between inner and outer radius centered at center
        must be odd elements and square
        if inner = 0 then no inner bound
    
        return values, (xinds, yinds) from original input array
        """    
        daRng = numpy.arange(0,data.shape[0],1)
        xGrid, yGrid = numpy.meshgrid(daRng, daRng)
        xOff = xGrid - center[0]
        yOff = yGrid - center[1]
        dist = numpy.sqrt(xOff**2+yOff**2)
        inds = numpy.nonzero((dist <= outer) & (dist >=inner))
        values = data[inds]
        return values, inds
        
    def getSkyCounts(self, denseData, denseCenter):
        """Determine the sky value using the median value determined in the outer annulus
        
        Returns:
        the median value within the outer annulus
        """

        values, inds = radialExtract(denseData, denseCenter, 
                                     self.skyAnnulus[1] * self.gridDense,
                                     self.skyAnulus[0] * self.gridDense)
        values
        
    def getPSFCounts(self, denseData, denseCenter):
        """Determine the sky value using the median value determined in the outer annulus
        
        Returns:
        the median value within the outer annulus
        """
        # Annulus was given in pixels corresponding to the original image data
        # since then we have extracted a region and densified it, so the Annulus
        # values must be scaled accordingly
        values, inds = radialExtract(denseData, denseCenter, 
                                     self.inrad * self.gridDense)
        values
        
        
class DiffPhotObj(object):
    """An object containing differential photometry information
    """
    def __init__(self, img, targCentroid, compCentroids, inrad, skyAnnulus, resolution, spline):
        """inputs:
        img: a baseImg or subclass of
        targCentroid: PyGuide centoid object for target star
        compCoords: list of PyGuide centroid objects for comparison stars
        """
        self.img = img
        self.targCentroid = targCentroid # PyGuide centriod
        self.compCentroids = compCentroids # list of PyGuide centroids
        apPhot = ApPhot(
            data = img.data, 
            inrad = inrad,
            skyAnnulus = skyAnnulus,
            gridDense = resolution,
            splineOrder = spline)
        # do photometry for target, offset pyguide to 0,0.
        if self.targCentroid.isOK:
            self.targPhot = apPhot.fuckinDoIt(numpy.asarray(self.targCentroid.xyCtr) - 0.5)
        else:
            self.targPhot = None
        self.compPhot = []
        for comp in self.compCentroids:
            if comp.isOK:
                self.compPhot.append(apPhot.fuckinDoIt(numpy.asarray(comp.xyCtr) - 0.5))
            else:
                self.compPhot.append(None)

    @property
    def targCounts(self):
        """photometry extracted target counts, normalized by exposure time
        Returns numpy.nan if no object was found by PyGuide at the specified location
        """
        if self.targPhot: # targPhot = None if pyguide centroid is not ok
            return self.targPhot.counts/self.img.exptime
        else:
            return numpy.nan
        
    @property
    def compCounts(self):
        """photometry extracted array of comparision counts, normalized by exposure time
        Returns numpy.nan if no object was found by PyGuide at the specified location
        """
        out = []
        for comp in self.compPhot:
            if comp: # comp=None if pyguide centroid is not ok.
                out.append(comp.counts/self.img.exptime)
            else: 
                out.append(numpy.nan)
        return out
        
    @property
    def diffCounts(self):
        """return the target counts normalized by the comparison counts summed
        """
        return self.targCounts / numpy.sum(self.compCounts)
                                     
                
        