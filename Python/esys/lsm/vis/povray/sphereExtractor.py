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
from .color        import RgbColor
from .sphere       import Sphere
import string

class SphereExtractor(core.SphereExtractor):
    def __init__(
        self,
        radiusMap   = lambda dataRecord: dataRecord.getRadius(),
        centerMap   = lambda dataRecord: dataRecord.getCenter(),
        modifierMap = lambda dataRecord: RgbColor(1,1,1),
        radiusScale = 1.0
    ):
        core.SphereExtractor.__init__(
            self,
            radiusMap,
            centerMap,
            modifierMap,
            radiusScale
        )

    def getGlyph(self, record):
        """
        Returns the sphere associated with the data in C{record}.
        """
        glyph = \
            Sphere(
                radius=self.getRadius(record),
                center=self.getCenter(record)
            )
        glyph.apply(self.getModifier(record))
        return glyph

    def writeSdl(self, f, record):
        self.getGlyph(record).writeSdl(f)




