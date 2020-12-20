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
Defines the L{SphereExtractor} class.
"""
from .color import Colors
from .extractor import Extractor

class SphereExtractor(Extractor):
    """
    Objects of this class can be used in conjunction with the
    L{esys.lsm.vis.core.GlyphData} class for extracting sphere
    info from data-records.
    """
    def __init__(
        self,
        radiusMap   = lambda dataRecord: dataRecord.getRadius(),
        centerMap   = lambda dataRecord: dataRecord.getCenter(),
        modifierMap = lambda dataRecord: Colors.Red,
        radiusScale = 1.0
    ):
        """
        Constructs the extractor.
        @type radiusMap: callable
        @param radiusMap: A callable which accepts a single data-record
        argument and returns a radius (float) value.
        @type centerMap: callable
        @param centerMap: A callable which accepts a single data-record
        argument and returns a 3 float-element sequence (ie a 3D coordinate).
        @type modifierMap: callable
        @param modifierMap: A callable which accepts a single data-record
        argument and returns an object modifier (or sequence of modifiers).
        @type radiusScale: float
        @param radiusScale: Scaling factor applied to all sphere radii.
        """
        self.radiusMap   = radiusMap
        self.centerMap   = centerMap
        self.modifierMap = modifierMap
        self.radiusScale = radiusScale

    def getRadiusScale(self):
        return self.radiusScale

    def getRadius(self, dataRecord):
        return self.radiusMap(dataRecord)*self.getRadiusScale()

    def getCenter(self, dataRecord):
        return self.centerMap(dataRecord)

    def getModifier(self, dataRecord):
        return self.modifierMap(dataRecord)
