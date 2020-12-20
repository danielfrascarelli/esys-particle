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

class Arrow(object):
    """
    Objects of this class represent arrows with conical head
    and cylindrical tail.
    """
    def __init__(self, tailPt, headPt, tailRadius, headRadius, headLength):
        """
        Initialises arrow.
        """
        self.headPt     = tuple(headPt)
        self.tailPt     = tuple(tailPt)
        self.tailRadius = tailRadius
        self.headRadius = headRadius
        self.headLength = headLength

    def getTailPt(self):
        """
        Returns the coordinate of the tail of this arrow.
        @return: Tail coordinate of this arrow.
        """
        return self.tailPt

    def getHeadPt(self):
        """
        Returns the coordinate of the pointy head of this arrow.
        @return: Head coordinate of this arrow.
        """
        return self.headPt

    def getHeadRadius(self):
        """
        Returns the max radius of the head of this arrow.
        @return: Radius of base of this arrow's conical head.
        """
        return self.headRadius

    def getTailRadius(self):
        """
        Returns the radius of the cylindrical tail of this arrow.
        @return: Radius of this arrow's cylindrical tail.
        """
        return self.tailRadius

    def getHeadLength(self):
        """
        Returns the height of the cone which forms the head of this arrow.
        @return: Height of cone which forms the head of this arrow.
        """
        return self.headLength

    def getLength(self):
        """
        Returns the distance between this arrows head-point and tail-point.
        @rtype: float
        @return: length of this arrow from head-point to tail-point.
        """
        return (Vec3(self.getHeadPt())-Vec3(self.getTailPt())).norm()



