"""Image objects and ways to combine them

todo?:
-add a reject algorithms to imcombine
-how to deal with an unfull calibrator? ie without any flats, just bias?
"""
import pyfits
import numpy

class Img(object):
    """Base class for images
    """
    def __init__(self, path, dateObs, exptime, filter=None, type=None, calibrator=None):
    """    
    inputs: 
    path: path to the image file
    dateObs: date of observation, a python datetime object
    exptime: the exposure time
    type: image type, defined by subclasses...
    filter: filter used
    """
        self.path = path
        self.dateObs = dateObjs
        self.exptime = float(exptime)
        if type not in ['bias', 'flat', 'light']:
            raise RuntimeError('image must be of the type: bias, flat, or light, got %s' % type)
        self.type = type
        self.filter = filter
        self.calibrator = calibrator
        
    @property        
    def data(self):
        img = pyfits.open(self.filename)
        data = img[0].data
        img.close()
        if self.calibrator:
            # apply calibrations
            data = self.calibrator.calibrate(data)
        return data
    
    def setCalibrator(self, calibrator):
        """Add (or switch) a calibrator, not sure if this will be useful
        """
        self.calibrator = calibrator
        
class Bias(Img):
    """Raw bias image
    """
    def __init__(self, path, dateObs, exptime=0, type='bias'):
        Img.__init__(self, path, dateObs, exptime, type)

class Flat(Img):
    """Raw Flat-field image
    """
    def __init__(self, path, dateObs, exptime, filter, type='flat'):
        Img.__init__(self, path, dateObs, exptime, filter, type)

class Light(Img):
    """Raw exposure (probably of astronomical objects)
    """
    def __init__(self, path, dateObs, exptime, filter, type='light'):
        Img.__init__(self, path, dateObs, exptime, filter, type)


class ImgCombine(object):
    """A generic method for combining images
    returns a combined image
    """
    def __init__(self, combType=None):
        """inputs:
        combType: a numpy.<something> function
        """
        self.combType = combType
    
    def __call__(imgList)
        """imgList: a list of Img, or subclass objects
        """
        nImg = len(imgList)
        x,y = imgList[0].shape
        # make a 3d image stack, then operate on the 3rd axis
        # initialize as zeros
        imgStack = numpy.zeros((x, y, nImg))
        for iter, img in enumerate(imgList):
            imgStack[:,:,iter]=img.data
        return self.combtype(imgStack, axis=2)

zeroCombine = ImgCombine(numpy.average)
flatCombine = ImgCombine(numpy.median)

#calibrator should calibrate the img

class Calibrator(object):
    """object that holds the calibration frames
    """
    def __init__(self, biasList, flatList):
        """
        biasList: a list of img.Bias image objects. 
        flatList: a list of img.Flat Objects.
        """
        
        self.bias = zeroCombine(biasList) if biasList else None # combine into master bias
        self.flat = dealWithFlats(flatList) if flatList else None# dict of master flats indexed by filter
            
    def applyBias(self, imgData):
        """inputs:
        img: an img.Light object
        filter: an image filter, if None, no flat fielding is applied
        returns:
        calibrated image array
        """
        # subtract the bias
        return imgData - self.bias
        
    def applyFlat(self, imgData, filter):
        """inputs:
        img: an img.Light object
        filter: an image filter, if None, no flat fielding is applied
        returns:
        calibrated image array
        """
        return imgData / self.flat[filter]
    
    def calibrate(self, imgData):
        """Apply bias and flat field to imgData
        """
        if self.bias:
            imgData = self.applyBias(imgData)
            if self.flat:
                imgData = self.applyFlat(imgData)
        return imgData
         
        
            
    def dealWithFlats(self, flatList):
        """sort and combine flatList into a flat dictionary, indexed by the filter
        """
        # first loop through and discover all filters
        flatDict = {}
        for flat in flatList:
            filter = flat.filter
            if filter not in flatDict.keys():
                # add the flat key to the dictionary, initialize with empty list
                flatDict[filter] = []
            # add the image to the key, with bias removed
            flatDict[filter].append(flat -  self.bias)
        # next combine lists of flat images into a master flat for each filter
        for flat, imgList in flatDict.iteritems():
            masterFlat = flatCombine(imgList)
            # overwrite image list in dictionary with master flat
            # and subtract the bias from it
            flatDict[flat] = masterFlat 
        return flatDict
            
def imgLister(fileList, type, instDict, calibrator=None):
    """function for building image lists from fileLists
    inputs:
    filelist: a list of strings defining each file
    type: 'bias', 'flat', or 'light'
    instDict: eg camera.flareCam.  Holds header solutions, etc
    calibrator: an image calibrator to go along with
    
    output:
    imgList, a list of image objects, ordered by observed date
    """
    imgList = []
    for file in fileList:
        # extract and save useful header data from each image
        imgFile = pyfits.open(file)
        if type = 'bias':
            filter = None
        else:
            filter = imgFile[0].header[instDict['filter']].strip() # get filter used, strip whitespace
        dateObs = instDict.parseDate(imgFile[0].header[instDict['dateObs']])
        exptime = imgFile[0].header[instDict['exptime']]
        imgList.append(
            Img(
                path = file, dateObs = dateObs, exptime = exptime, 
                filter = filter, type = type, 
                calibrator = calibrator,
            )
        )
        imgFile.close()           
        
           
        
        
        
        
        
    
