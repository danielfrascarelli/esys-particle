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

"""Defines the L{TruncatedCone} and L{Cone} base classes"""
from .exception import raiseNotImplemented
from esys.lsm.util import Vec3

class TruncatedCone(object):
    """
    Objects of this class represent truncated cones.
    """
    def __init__(self, endPt1, endPt2, radius1, radius2=0):
        """
        Initialises cone with center-end-point coordinates and
        corresponding radii.
        @type endPt1: sequence of 3 floats
        @param endPt1: Center coordinate of one end.
        @type endPt2: sequence of 3 floats
        @param endPt2: Center coordinate of other end.
        @type radius1: float
        @param radius1: End 1 radius.
        @type radius2: float
        @param radius2: End 2 radius.
        """
        self.endPt1 = tuple(endPt1)
        self.endPt2 = tuple(endPt2)
        self.radius1 = radius1
        self.radius2 = radius2

    def getEndPt1(self):
        """
        Returns the coordinate of the center coordinate of end 1.
        @return: Center coordinate of end 1.
        """
        return self.endPt1

    def getEndPt2(self):
        """
        Returns the coordinate of the center coordinate of end 2.
        @return: Center coordinate of end 2.
        """
        return self.endPt2

    def getCenter(self):
        """
        Returns the coordinate of the center of this cone.
        @return: Center coordinate of this cone.
        """
        return (Vec3(self.getEndPt1()) + Vec3(self.getEndPt2()))*0.5

    def getHeight(self):
        """
        Returns the height of this cone.
        @rtype: float
        @return: Center coordinate of this cone.
        """
        return (Vec3(self.getEndPt1()) - Vec3(self.getEndPt2())).norm()

    def getRadius1(self):
        """
        Returns the end 1 radius of this cone.
        @return: End 1 radius.
        """
        return self.radius1

    def getRadius2(self):
        """
        Returns the end 2 radius of this cone.
        @return: End 2 radius.
        """
        return self.radius2

class Cone(TruncatedCone):
    """
    Objects of this class represent cones.
    """
    def __init__(self, basePt, tipPt, radius):
        """
        Initialises cone with center-end-point coordinates and
        corresponding radii.
        @type basePt: sequence of 3 floats
        @param basePt: Center coordinate of one end.
        @type tipPt: sequence of 3 floats
        @param tipPt: Center coordinate of other end.
        @type radius: float
        @param radius: Radius at C{basePt}.
        """
        TruncatedCone.__init__(self, basePt, tipPt, radius, 0)

    def getBasePt(self):
        return self.getEndPt1()

    def getTipPt(self):
        return self.getEndPt2()

    def getRadius(self):
        return self.getRadius1()






