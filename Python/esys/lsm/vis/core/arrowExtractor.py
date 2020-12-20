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

"""Defines the L{Arrow} base class"""
from .exception import raiseNotImplemented
from esys.lsm.util import Vec3

class ArrowExtractor(object):
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
        self.tailPtMap     = tailPtMap
        self.vecMap        = vecMap
        self.tailRadiusMap = tailRadiusMap
        self.headRadiusMap = headRadiusMap
        self.headLengthMap = headLengthMap
        self.modifierMap   = modifierMap
        self.lengthScale   = lengthScale

    def getTailPt(self, record):
        return self.tailPtMap(record)

    def getVec(self, record):
        return self.vecMap(record)

    def getHeadPt(self, record):
        return \
            (
                Vec3(self.getTailPt(record))
                +
                (Vec3(self.getVec(record))*self.getLengthScale())
            )

    def getHeadRadius(self, record):
        return self.headRadiusMap(record)

    def getTailRadius(self, record):
        return self.tailRadiusMap(record)

    def getHeadLength(self, record):
        return self.headLengthMap(record)

    def getModifier(self, record):
        return self.modifierMap(record)

    def getLengthScale(self):
        return self.lengthScale

