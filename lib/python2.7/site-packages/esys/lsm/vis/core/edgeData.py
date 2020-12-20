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
Defines the L{EdgeData} class.
"""
from .edgeExtractor import EdgeExtractor

class EdgeData(object):
    """
    Used to render data which contains piecewise linear I{edges}.
    """
    def __init__(self, data, edgeExtractor):
        """
        Constructor.
        @type data: iterable
        @param data: Record-data elements containing piecewise-linear
        coordinates.
        @type edgeExtractor: L{EdgeExtractor}
        @param edgeExtractor: The object used to extract edge-endpoints
        from the C{data}.
        """
        self.data          = data
        self.edgeExtractor = edgeExtractor

    def getData(self):
        """
        Returns the I{iterable} data associated with this C{EdgeData} object.
        @rtype: object
        @return: The I{iterable} object which contains the end-points
        of piecewise linear edges.
        """
        return self.data

    def getEdgeExtractor(self):
        """
        Returns the I{extractor} object used to obtain end-points from
        the data-records contained in the C{self.getData()} iterable.
        @rtype: L{EdgeExtractor}
        @return: L{EdgeExtractor} object used to obtain edge-points data points.
        """
        return self.edgeExtractor
