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

"""Defines the L{Cylinder} and L{Disk} base classes"""
from __future__ import division
from .exception import raiseNotImplemented
from esys.lsm.util import Vec3

class Cylinder(object):
    """
    Objects of this class represent cylinders.
    """
    def __init__(self, endPt1, endPt2, radius):
        """
        Initialises cylinder with center-end-point coordinates and radius.
        @type endPt1: sequence of 3 floats
        @param endPt1: Center coordinate of one end.
        @type endPt2: sequence of 3 floats
        @param endPt2: Center coordinate of other end.
        @type radius: float
        @param radius: Radius of cylinder.
        """
        self.endPt1 = tuple(endPt1)
        self.endPt2 = tuple(endPt2)
        self.radius = radius

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
        Returns the coordinate of the center of this cylinder.
        @return: Center coordinate of this cylinder.
        """
        return (Vec3(self.getEndPt1()) + Vec3(self.getEndPt2()))*0.5

    def getHeight(self):
        """
        Returns the height of this cylinder.
        @rtype: float
        @return: Center coordinate of this cylinder.
        """
        return (Vec3(self.getEndPt1()) - Vec3(self.getEndPt2())).norm()

    def getRadius(self):
        """
        Returns the radius of this cylinder.
        @return: Radius of this cylinder.
        """
        return self.radius

class Disk(Cylinder):
    """
    Objects of this class represent cylinders.
    """
    def __init__(self, center, radius, height, direction=None):
        if (direction == None):
            direction = core.Vec3(0,0,1)
        direction = core.Vec3(direction)
        nrm = direction.norm()
        if (nrm == 0.0):
            direction = core.Vec3(0,0,1)
            nrm = 1.0

        center = core.Vec3(centre)
        direction = direction/nrm
        Cylinder.__init__(
            self,
            endPt1 = center - (direction*(height*0.5)),
            endPt2 = center + (direction*(height*0.5)),
            radius = radius
        )
