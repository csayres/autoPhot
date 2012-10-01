"""For generating triangles and determining offsets
"""

import numpy
import PyGuide
import numpy.linalg
import itertools


class TriangleHash(object):
    """An object representing a bunch of triangles
    """
    def __init__(self, coordList):
        """Inputs:
        - coordList     - a 2d numpy array xy positions of found stars 
            on ccd (output from PyGuide.findStars)
        """
        self.coordList = coordList
        self.verts = None
        self.triSpace = None

        
    def arrangeVerts(self, verts):
        """Inputs:
        verts - a 3x2 numpy array containing un-arranged triangle vertices
        
        Output:
        a 3x2 numpy array with vertices arranged as so:
            A and B define the
            shortest side, B and C the longest side, and A and C define
            the intermediate length side
        """
        # moves last element to front, shifts rest
        vertsShifted = numpy.roll(verts, 1, axis=0) 
        # compute distances
        dist = numpy.sqrt(numpy.sum( (verts - vertsShifted)**2 ,1) )
        # match up between shifted and unshifted vertices, the way distances were determined
        indMap = numpy.array([ [0, 2], [1, 0], [2, 1] ]) 
        distSortInd = numpy.argsort(dist) # returns indices to sort by
        sortedPairs = indMap[distSortInd, :] # sorted pairs, increasing distance
        
        # find vertex between short and medium length side: A
        oneOf = sortedPairs[0:2, :] # shortest 2 lengths.  Common index will be A
        uniqueInds, allInds = numpy.unique(oneOf, return_inverse = True)
        # whats going on here?    
#         >>> g
#         array([[0, 1],
#                [1, 2]])
#         >>> numpy.unique(g, return_inverse=True)
#         (array([0, 1, 2]), array([0, 1, 1, 2]))
#       the difference in the sum of both of these arrays is 1 (the repeated number)
        AInd = numpy.sum(allInds) - numpy.sum(uniqueInds)
        A = verts[AInd]
        
        # B vertex will connect with A to make the short side
        shortSide = sortedPairs[0, :] # shortest side indices, contains AInd
        BInd = numpy.sum(shortSide) - AInd # sum and subtract AInd to reveal BInd
        B = verts[BInd]
        
        # C vertex is the only remaining
        # possible Inds = [0, 1, 2], sum is three
        CInd = 3 - (AInd + BInd)
        C = verts[CInd]
       
        return numpy.vstack((A, B, C))
        
    def triangleMatch(self, triSpace1, triSpace2):
        """Match two sets of lists within a 2% tolerance that have been 
        converted to triangle space arrays
        Inputs:
        triSpace1: output from triangleSpace() n x 2 array with columns x_t, y_t
        triSpace2: output from triangleSpace() n x 2 array with columns x_t, y_t
        
        Output:
        n x 2 array containing indices of potential matches between
        triSpace1 and triSpace2.
        
        note: more than 1 object may be matched to an object in triSpace2
        and vice versa.  triangleVerify() will ultimitely decide the best
        matches in real space.
        """        
        tol = 0.02 # must match to 2 percent
        potentialMatches = []
        # sort triangles by their product of x_t and y_t (most unique first)
        triSort1 = numpy.argsort(triSpace1[:,0] * triSpace1[:,1])
        triSort2 = numpy.argsort(triSpace2[:,0] * triSpace2[:,1])
        
        for ind1 in triSort1:
            for ind2 in triSort2:
                ratio = numpy.abs(1 - (triSpace1[ind1, :] / triSpace2[ind2, :]))
                if numpy.max(ratio) < tol:
                    potentialMatches.append([ind1, ind2])
                if len(potentialMatches) > 10:
                    break # if we have more than 10 potential triangles
                    
            # determine % variation between y_t values
#             yratio = numpy.abs(1 - (tri1[1] / triSpace2[:, 1])) # abs because tol is +-2%
#             tolInd = numpy.nonzero(yratio < tol) # keep only those which y_t tol < .02
#             for ind2 in tolInd:
#                 # now x_t must also be within 2%
#                 print 'tri1.shape: ', tri1.shape, 'tri
#                 xratio = numpy.abs(1 - (tri1[0] / triSpace2[ind2, 0]))
#                 print 'xratio: ', xratio
#                 if xratio < tol:
#                     potentialMatches.append([ind1, ind2])
                    
                    
                    
        return numpy.asarray(potentialMatches)
            
    def triangleVerify(self, potentialMatches, triVert1, triVert2):
        """Take potential matches from the triSpace matching done in triangleMatch()
        and check to see if they match up in real space!
        
        Inputs: 
        potentialMatches - output of triangleMatch
        triVert1 - output of makeTriVerts() 3D array n x 3pts x 2coords
        triVert2 - output of makeTriVerts() 3D array n x 3pts x 2coords
        
        Output:
        n x 2 array with match indices between triVert1 and triVert2        
        """
        tol = 2 # pixels
        tri1Matched = [] # holds indices already matched, can't be reused
        tri2Matched = []
        matched = []
        for match in potentialMatches:
            # check to see if either match is already assigned
            if match[0] in tri1Matched or match[1] in tri2Matched:
                continue
            # normalize by vertex A
            tri1 = triVert1[match[0], :, :] #2D
            tri1 = tri1 - tri1[0, :]
            tri2 = triVert2[match[1], :, :]
            tri2 = tri2 - tri2[0, :]
            diff = tri1 - tri2 # first vertex at 0,0 due to normalization
            dist = numpy.sqrt(numpy.sum(diff**2, 1)) # distance between verts
            if numpy.max(dist) < tol:
                # both unnormalized vertices are less than our pixel tolerance
                tri1Matched.append(match[0])
                tri2Matched.append(match[1])
                matched.append([match[0], match[1]])
        #print 'returned best val: ', numpy.hstack((tri1Matched, tri2Matched))
        return numpy.asarray(matched)
        
    def getOffset(self, verifiedMatches, triVert1, triVert2):
        """Take verified matches from the triangleVerify() matching
        
        Inputs: 
        verifiedMatches - output of triangleVerify()
        triVert1 - output of makeTriVerts() 3D array n x 3pts x 2coords
        triVert2 - output of makeTriVerts() 3D array n x 3pts x 2coords
        
        Output:
        2 element array determining the offset to be applied       
        """             
        # just take the average offset of all for now
        offsets = []
        for match in verifiedMatches:
            tri1 = triVert1[match[0], :, :]
            tri2 = triVert2[match[1], :, :]
            diff = numpy.subtract(tri1, tri2)
            offsets.append(numpy.mean(diff, 0))
        return numpy.mean(numpy.asarray(offsets), 0)
                                  
    def makeTriVerts(self):
        """Return a list of triangle Vertices, arranged by side length
        see: arrangeVerts()
        
        Output: 3D array
        """
        triangles = map(self.arrangeVerts, itertools.combinations(self.coordList, 3))
#         triangles = []
#         old way, slightly slower:                            
#         for ind1, coord1 in enumerate(self.coordList):
#             for ind2, coord2 in enumerate(self.coordList):
#                 if ind2 <= ind1: continue
#                 for ind3, coord3 in enumerate(self.coordList):
#                     if ind3 <= ind2: continue
#                     ABC = self.arrangeVerts(numpy.array([coord1, coord2, coord3]))
#                     triangles.append(ABC)

        triangles = numpy.asarray(triangles) # 3D array
        return triangles        

    def triangleSpace(self, triangles):
        """Create a output triangle space to search through
        
        Inputs: 
        triangles - A 3D numpy array as is output from makeTriVerts (ABC verts in a specific order)
            Each element (triangle), contains 3 vertices, containing an x and y point:
            nx3x2 array where n is num of triangles
            
        Outputs:
        triangleSpace - A 2D numpy array n x (x_t, y_t)
            x_t = CB dot CA
            y_t = len(A) / len(C)
        """
        out = []
        for tri in triangles:
            CA = tri[2] - tri[0]
            CB = tri[2] - tri[1] # longest side
            BA = tri[1] - tri[0] # shortest side
            x_t = numpy.dot(CB, CA)
            y_t = numpy.linalg.norm(CB) / numpy.linalg.norm(BA)
            out.append(numpy.array([x_t, y_t]))
        return numpy.asarray(out)

    def setup(self):
        # hack for now to keep computations out of outer loop
        # unless they are needed
        self.verts = self.makeTriVerts()
        self.triSpace = self.triangleSpace(self.verts)
    
    # could move this work to a standalone object
    def hashItOut(self, otherHash):
        """Take the current hash compare it with another, and return the 
        determined offset
        """
        print 'hashin it out'
        self.setup()
        otherHash.setup() # ugh
        potentialMatches = self.triangleMatch(self.triSpace, otherHash.triSpace)
        verifiedMatches = self.triangleVerify(potentialMatches, self.verts, otherHash.verts)
        return self.getOffset(verifiedMatches, self.verts, otherHash.verts)
            