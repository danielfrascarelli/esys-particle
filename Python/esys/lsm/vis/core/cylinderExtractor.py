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
Defines the L{CylinderExtractor} class.
"""
from .color import Colors
from .extractor import Extractor

class CylinderExtractor(Extractor):
    """
    Objects of this class can be used in conjunction with the
    L{esys.lsm.vis.core.GlyphData} class for extracting cylinder
    info from data-records.
    """
    def __init__(
        self,
        radiusMap   = lambda dataRecord: dataRecord.getRadius(),
        endPt1Map   = lambda dataRecord: dataRecord.getEndPt1(),
        endPt2Map   = lambda dataRecord: dataRecord.getEndPt2(),
        modifierMap = lambda dataRecord: None,
        radiusScale = 1.0
    ):
        """
        Constructs the extractor.
        @type radiusMap: callable
        @param radiusMap: A callable which accepts a single data-record
        argument and returns a radius (float) value.
        @type endPt1Map: callable
        @param endPt1Map: A callable which accepts a single data-record
        argument and returns a 3 float-element sequence (ie a 3D coordinate).
        @type endPt2Map: callable
        @param endPt2Map: A callable which accepts a single data-record
        argument and returns a 3 float-element sequence (ie a 3D coordinate).
        @type modifierMap: callable
        @param modifierMap: A callable which accepts a single data-record
        argument and returns an object modifier (or sequence of modifiers).
        """
        self.radiusMap   = radiusMap
        self.endPt1Map   = endPt1Map
        self.endPt2Map   = endPt2Map
        self.modifierMap = modifierMap
        self.radiusScale = radiusScale

    def getRadius(self, dataRecord):
        return self.radiusMap(dataRecord)*self.radiusScale

    def getEndPt1(self, dataRecord):
        return self.endPt1Map(dataRecord)

    def getEndPt2(self, dataRecord):
        return self.endPt2Map(dataRecord)

    def getModifier(self, dataRecord):
        return self.modifierMap(dataRecord)

    def getRadiusScale(self):
        return self.radiusScale

class DiskExtractor(Extractor):
    """
    Objects of this class can be used in conjunction with the
    L{esys.lsm.vis.core.GlyphData} class for extracting cylinder
    info from data-records.
    """
    def __init__(
        self,
        radiusMap    = lambda dataRecord: dataRecord.getRadius(),
        centerMap    = lambda dataRecord: dataRecord.getCenter(),
        heightMap    = lambda dataRecord: 0.001,
        directionMap = lambda dataRecord: None,
        modifierMap  = lambda dataRecord: None,
        radiusScale  = 1.0
    ):
        """
        Constructs the extractor.
        @type radiusMap: callable
        @param radiusMap: A callable which accepts a single data-record
        argument and returns a radius (float) value.
        @type centerMap: callable
        @param centerMap: A callable which accepts a single data-record
        argument and returns a 3 float-element sequence (ie a 3D coordinate).
        @type heightMap: callable
        @param heightMap: A callable which accepts a single data-record
        argument and returns a float (height of the cylinder).
        @type modifierMap: callable
        @param modifierMap: A callable which accepts a single data-record
        argument and returns a color.
        """
        self.radiusMap    = radiusMap
        self.centerMap    = centerMap
        self.heightMap    = heightMap
        self.directionMap = directionMap
        self.modifierMap  = modifierMap
        self.radiusScale  = radiusScale

    def getRadius(self, dataRecord):
        return self.radiusMap(dataRecord)*self.radiusScale

    def getCenter(self, dataRecord):
        return self.centerMap(dataRecord)

    def getHeight(self, dataRecord):
        return self.heightMap(dataRecord)

    def getDirection(self, dataRecord):
        return self.directionMap(dataRecord)

    def getModifier(self, dataRecord):
        return self.modifierMap(dataRecord)

    def getRadiusScale(self):
        return self.radiusScale

