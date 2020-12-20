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
from .sphere       import Sphere

class SphereExtractor(core.SphereExtractor):
    def __init__(
        self,
        radiusMap   = lambda dataRecord: dataRecord.getRadius(),
        centerMap   = lambda dataRecord: dataRecord.getCenter(),
        modifierMap = lambda dataRecord: None,
        radiusScale = 1.0
    ):
        core.SphereExtractor.__init__(
            self,
            radiusMap,
            centerMap,
            modifierMap,
            radiusScale
        )

    def addActor(self, record, addTo):
        """
        Adds the sphere-actor associated with the data in
        C{record} to C{addTo}.
        """
        glyph = \
            Sphere(
                center=self.getCenter(record),
                radius=self.getRadius(record)
            )
        glyph.apply(self.getModifier(record))
        glyph.addActor(addTo)

