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
from .arrow        import Arrow

class ArrowExtractor(core.ArrowExtractor):
    def __init__(
        self,
        tailPtMap,
        vecMap,
        tailRadiusMap = lambda record: 0.10,
        headRadiusMap = lambda record: 0.25,
        headLengthMap = lambda record: 0.75,
        modifierMap   = lambda record: None,
        lengthScale   = 1.0
    ):
        core.ArrowExtractor.__init__(
            self,
            tailPtMap,
            vecMap,
            tailRadiusMap,
            headRadiusMap,
            headLengthMap,
            modifierMap,
            lengthScale
        )

    def addActor(self, record, addTo):
        """
        Adds the arrow-actor associated with the data in
        C{record} to C{addTo}.
        """
        glyph = \
            Arrow(
                tailPt     = self.getTailPt(record),
                headPt     = self.getHeadPt(record),
                tailRadius = self.getTailRadius(record),
                headRadius = self.getHeadRadius(record),
                headLength = self.getHeadLength(record)
            )
        glyph.apply(self.getModifier(record))
        glyph.addActor(addTo)
