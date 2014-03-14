"""Output Results in various ways
"""
import viz
import matplotlib.pyplot as plt
import itertools
import numpy

ZEROPOINT = 20

class Dump(object):
    """A class for dumping information to a directory
    """
    def __init__(self, outDir, driver, units = 'Flux'):
        """
        Inputs:
        outDir      an output directory
        driver      a Driver object
        units       'Flux' or 'Mag'
        """
        self.outDir = outDir
        self.driver = driver
        if units not in ['Flux', 'Mag']:
            raise RuntimeError('units must be "Flux" or "Mag", received: %s' % units)
        self.units = units
        
    def saveField(self):
        """generate and save a figure from the extraction field
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # note, field solution is the last extracted image
        viz.showFieldSolution(ax, self.driver.fieldSolution)
        figName = 'field.png'
        fig.savefig(self.outDir + figName)
        
    def saveDatFile(self):
        targCounts, compCounts, diffCounts = self.getPhotData()
        nComps = compCounts.shape[1]
        fileLines = [
            'Aperture (radius, pixels) = %i' % self.driver.config.phot.inrad,
            'Sky Annulus (inner and outer radius, pixels) = %i, %i' % (self.driver.config.phot.skyann[0],self.driver.config.phot.skyann[1]),
            '------------------------------',
            ] 
        lowerHead = ['Normalized Differential Photometry (Units = %s), Target (Units = %s), ' % (self.units, self.units)]
        for num in range(nComps):
            lowerHead.append('Comparison%i (Units = %s), ' % (num, self.units))
        lowerHead = ''.join(lowerHead)
        fileLines.append(lowerHead)
        for diff,targ,comp in itertools.izip(diffCounts, targCounts, compCounts):
            slam = [diff]
            slam.append(targ)
            slam.extend(comp)
            dataLine = []
            for phot in slam:
                dataLine.append('%f' % phot)
            dataLine = ', '.join(dataLine)
            fileLines.append(dataLine)
        
        # now write the file
        fileName = 'photOut.dat'
        with open(self.outDir + fileName, 'w') as f:
            for line in fileLines:
                f.write(line+'\n')
    
    def saveLightCurves(self):
        """print out all light curves
        """
        targCounts, compCounts, diffCounts = self.getPhotData()
        fig = plt.figure()
        
    
    def getPhotData(self):
        """Loop through each extraction get all the fwhms
        """
        targCounts = [] #1D
        compCounts = [] #2D
        diffCounts = [] #1D       
        
        for photObj in self.driver.df:
            targCounts.append(photObj.targCounts)
            compCounts.append(photObj.compCounts)  
            diffCounts.append(photObj.diffCounts)
            
        targCounts = numpy.asarray(targCounts) #1D
        compCounts = numpy.asarray(compCounts) #2D
        diffCounts = numpy.asarray(diffCounts) #1D
        
        if self.units == 'Mag':
            # convert counts to and (arbitrary magnitude)
            targCounts = self.flux2mag(targCounts)
            compCounts = self.flux2mag(compCounts)
            diffCounts = self.flux2mag(diffCounts)
        
        return targCounts, compCounts, diffCounts
    
    def flux2mag(self, array, zeropoint = ZEROPOINT):
        """convert flux to magnitude zeropoint is arbitrary
        """
        return -2.5*numpy.log10(array) + zeropoint        
            