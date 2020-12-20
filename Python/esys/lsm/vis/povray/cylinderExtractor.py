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

from esys.lsm.vis import core
from .cylinder     import Cylinder, Disk

import string

class CylinderExtractor(core.CylinderExtractor):
    def __init__(
        self,
        radiusMap   = lambda dataRecord: dataRecord.getRadius(),
        endPt1Map   = lambda dataRecord: dataRecord.getEndPt1(),
        endPt2Map   = lambda dataRecord: dataRecord.getEndPt2(),
        modifierMap = lambda dataRecord: None,
        radiusScale = 1.0
    ):
        core.CylinderExtractor.__init__(
            self,
            radiusMap,
            endPt1Map,
            endPt2Map,
            modifierMap,
            radiusScale
        )

    def getGlyph(self, record):
        """
        Returns the cylinder associated with the data in C{record}.
        """
        glyph = \
            Cylinder(
                endPt1=self.getEndPt1(record),
                endPt2=self.getEndPt2(record),
                radius=self.getRadius(record)
            )
        glyph.apply(self.getModifier(record))
        return glyph

    def writeSdl(self, f, record):
        self.getGlyph(record).writeSdl(f)

class DiskExtractor(core.DiskExtractor):
    def __init__(
        self,
        radiusMap    = lambda dataRecord: dataRecord.getRadius(),
        centerMap    = lambda dataRecord: dataRecord.getCenter(),
        heightMap    = lambda dataRecord: 0.001,
        directionMap = lambda dataRecord: None,
        modifierMap  = lambda dataRecord: None,
        radiusScale  = 1.0
    ):
        core.DiskExtractor.__init__(
            self,
            radiusMap,
            centerMap,
            heightMap,
            directionMap,
            modifierMap,
            radiusScale
        )

    def getGlyph(self, record):
        """
        Returns the cylinder associated with the data in C{record}.
        """
        glyph = \
            Disk(
                center    = self.getCenter(record),
                radius    = self.getRadius(record),
                height    = self.getHeight(record),
                direction = self.getDirection(record)
            )
        glyph.apply(self.getModifier(record))
        return glyph

    def writeSdl(self, f, record):
        self.getGlyph(record).writeSdl(f)

