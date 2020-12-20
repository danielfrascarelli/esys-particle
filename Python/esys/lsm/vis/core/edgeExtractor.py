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
Defines the L{EdgeExtractor} class.
"""
from .extractor import Extractor

class EdgeExtractor(Extractor):
    """
    Objects can be used in conjunction with L{EdgeData} objects
    to extract piecewise-linear edges from data-records. Each linear
    segment of an edge is rendered as a truncated cone, hence there
    is a radius associated with each point on the edge and a color
    associated with each linear segment on the edge.
    """
    def __init__(
        self,
        pointListMap    = lambda dataRecord: dataRecord.getPointList(),
        radiusListMap   = lambda dataRecord: 1.0,
        colorValListMap = lambda dataRecord: 0,
        radiusScale     = 1.0
    ):
        """
        Constructor.
        @type pointListMap: callable
        @param pointListMap: Callable which takes a single data-record
        as an argument and returns a list of coordinates representing
        the end-points of the linear segments of the edge.
        @type radiusListMap: callable
        @param radiusListMap: Callable which takes a single data-record
        as an argument and returns a list of radii values.
        @type colorValListMap: callable
        @param colorValListMap: Callable which takes a single data-record
        as an argument and returns a list of color values. Each color
        value is associated with a linear segment of the edge.
        @type radiusScale: float
        @param radiusScale: Scale which is multiplied by all radii.
        """
        self.radiusListMap   = radiusListMap
        self.pointListMap    = pointListMap
        self.colorValListMap = colorValListMap
        self.radiusScale     = radiusScale

    def getRadiusScale(self):
        """
        Returns the radius scaling factor.
        @rtype: float
        """
        return self.radiusScale

    def getRadiusList(self, dataRecord):
        """
        Returns the radii associated with each segment end-point.
        @rtype: sequence of floats
        @return: Sequence of float radii values.
        """
        return self.radiusListMap(dataRecord)

    def getPointList(self, dataRecord):
        """
        Returns the coordinates (3 float-element sequence)
        of each segment end-point.
        @rtype: sequence of coordinates
        @return: Sequence of coordinate values.
        """
        return self.pointListMap(dataRecord)

    def getColorValList(self, dataRecord):
        """
        Returns the color of each segment.
        @rtype: sequence of colors
        @return: Sequence of color objects.
        """
        return self.colorValListMap(dataRecord)
