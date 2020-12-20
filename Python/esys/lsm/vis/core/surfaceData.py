#############################################################
##                                                         ##
## Copyright (c) 2003-2017 by The University of Queensland ##
## Centre for Geoscience Computing                         ##
## http://earth.uq.edu.au/centre-geoscience-computing      ##
##                                                         ##
## Primary Business: Brisbane, Queensland, Australia       ##
## Licensed under the Open Software License version 3.0    ##
## http://www.apache.org/licenses/LICENSE-2.0              ##
##                                                         ##
#############################################################

"""
Defines the L{SurfaceData} class.
"""
from .pointExtractor import PointExtractor

class SurfaceData(object):
    """
    Represents a topographic surface.
    """
    def __init__(self, data, pointExtractor):
        """
        Constructs the object.
        @type data: iterable
        @param data: Object containing data which specifies discrete points
        on a surface.
        @type pointExtractor: L{PointExtractor}
        @param pointExtractor: An object to extract surface coordinates
          from the specified C{data}.
        """
        self.data           = data
        self.pointExtractor = pointExtractor

    def getData(self):
        """
        Returns the I{iterable} data associated with this C{SurfaceData} object.
        @rtype: object
        @return: The I{iterable} object which contains the surface coordinates.
        """
        return self.data

    def getPointExtractor(self):
        """
        Returns the I{extractor} object used to obtain points from
        the data-records contained in the C{self.getData()} iterable.
        @rtype: L{PointExtractor}
        @return: L{PointExtractor} object used to obtain surface data points.
        """
        return self.pointExtractor
