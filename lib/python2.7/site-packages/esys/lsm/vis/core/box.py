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

"""Defines the L{Box} base class"""

from esys.lsm.util import Vec3

class Box:
    """
    A rectangular box.
    """
    def __init__(self, minPt, maxPt):
        """
        Constructs axis-aligned box.
        @type minPt: sequence of 3 floats
        @param minPt: lower left back corner coordinate.
        @type maxPt: sequence of 3 floats
        @param maxPt: upper right front corner coordinate.
        """
        self.minPt = minPt
        self.maxPt = maxPt

    def getMinPt(self):
        return self.minPt

    def getMaxPt(self):
        return self.maxPt

    def getCenter(self):
        return (Vec3(self.getMaxPt())+Vec3(self.getMinPt()))*0.5

    def getSideLength(self):
        return (Vec3(self.getMaxPt())-Vec3(self.getMinPt()))
