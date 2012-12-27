"""Do differential photometry

note: PyGuide introduces a 0.5 pixel shift on centroiding (because pixels are
defined with centers at n*(1-0.5) rather than n

photometry set to nan if pyguide centroid is not ok--don't think it will ever
be not ok, but whatever.
"""
import phot
import numpy

class DiffPhotObj(object):
    """An object containing differential photometry information
    """
    def __init__(self, img, targCentriod, compCentriods, inrad, skyAnnulus, resolution, spline):
        """inputs:
        img: a baseImg or subclass of
        targCentroid: PyGuide centoid object for target star
        compCoords: list of PyGuide centroid objects for comparison stars
        """
        self.img = img
        self.targCentroid = targCentroid # PyGuide centriod
        self.compCentroids = compCentroids # list of PyGuide centroids
        apPhot = phot.ApPhot(
            data = img.data, 
            inrad = inrad,
            skyAnnulus = skyAnnulus,
            gridDense = resolution,
            splineOrder = spline)
        # do photometry for target, offset pyguide to 0,0.
        if self.targPhot.isOK:
            self.targPhot = apPhot.fuckinDoit(self.targCentroid.xyCtr - 0.5)
        else:
            self.targPhot = None
        self.compPhot = []
        for comp in self.compCentroids:
            if comp.isOK:
                self.compPhot.append(apPhot.fuckinDoit(comp.xyCtr - 0.5))
            else:
                self.compPhot.append(None)

    @property
    targCounts(self):
        """photometry extracted target counts, normalized by exposure time
        Returns numpy.nan if no object was found by PyGuide at the specified location
        """
        if self.targPhot: # targPhot = None if pyguide centroid is not ok
            return self.targPhot.counts/self.img.exptime
        else:
            return numpy.nan
        
    @property
    compCounts(self):
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
    diffCounts(self):
        """return the target counts normalized by the comparison counts summed
        """
        return self.targCounts / numpy.sum(self.compCounts)
        
        
        
        
    
            
            
                