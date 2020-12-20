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

"""Defines the L{Capsule} base class"""
from __future__ import division
from .exception import raiseNotImplemented
from esys.lsm.util import Vec3

class Capsule(object):
    """
    Objects of this class represent capsules, a cylinder capped with
    a hemisphere at each end.
    """
    def __init__(self, radius, endPt1=None, endPt2=None, segEndPt1=None, segEndPt2=None):
        """
        Initialises capsule with center-end-point coordinates and radius.
        @type radius: float
        @param radius: Radius of capsule.
        @type endPt1: sequence of 3 floats
        @param endPt1: Coordinate of one end ("apex" of hemisphere).
        @type endPt2: sequence of 3 floats
        @param endPt2: Center coordinate of other end ("apex" of hemisphere).
        @type segEndPt1: sequence of 3 floats
        @param segEndPt1: Coordinate of one end ("apex" of hemisphere).
        @type segEndPt2: sequence of 3 floats
        @param segEndPt2: Center coordinate of other end ("apex" of hemisphere).
        """
        self.radius = radius
        if (endPt1 != None) and (endPt2 != None):
            endPt1 = Vec3(endPt1)
            endPt2 = Vec3(endPt2)
            centre = 0.5*(endPt1 + endPt2)
            direction = (endPt2 - endPt1)
            length = direction.norm()
            direction /= length
            self.segEndPt1 = centre - (0.5*length)*direction
            self.segEndPt2 = centre + (0.5*length)*direction
        elif (segEndPt1 != None) and (segEndPt2 != None):
            self.segEndPt1 = Vec3(segEndPt1)
            self.segEndPt2 = Vec3(segEndPt2)
        else:
            raise Exception("Must specified end-points or segment-end-points, not both.")

    def getSegmentEndPt1(self):
        """
        Returns the coordinate of the center coordinate of end 1.
        @return: Center coordinate of end 1.
        """
        return self.segEndPt1

    def getSegmentEndPt2(self):
        """
        Returns the coordinate of the center coordinate of end 2.
        @return: Center coordinate of end 2.
        """
        return self.segEndPt2

    def getCenter(self):
        """
        Returns the coordinate of the center of this capsule.
        @return: Center coordinate of this capsule.
        """
        return (self.getSegmentEndPt1() + self.getSegmentEndPt2())*0.5

    def getSegmentDirection(self):
        """
        Returns the direction vector of the central axis.
        @return: Direction of longitudinal axis.
        """
        d = (self.getSegmentEndPt2() - self.getSegmentEndPt1())
        d /= d.norm()
        return d;

    def getLength(self):
        """
        Returns the length of this capsule.
        @rtype: float
        @return: Distance between hemi-sphere apexes.
        """
        return \
            2.0*self.radius() \
            +\
            (Vec3(self.getSegmentEndPt1()) - Vec3(self.getSegmentEndPt2())).norm()

    def getSegmentLength(self):
        """
        Returns the segment-length of this capsule.
        @rtype: float
        @return: minimum distance between end-points of cylinder.
        """
        return \
            (Vec3(self.getSegmentEndPt1()) - Vec3(self.getSegmentEndPt2())).norm()

    def getRadius(self):
        """
        Returns the radius of this capsule.
        @return: Radius of this capsule.
        """
        return self.radius

class CapsuleDisk:
    """
    Objects of this class represent 2-dimensional capsules.
    """
    def __init__(
        self,
        center,
        length,
        radius,
        height = None
    ):
        self.center = tuple(center)
        self.length = length
        self.radius = radius
        if (height == None):
            height = 0.05*radius
        self.height = height

    def getCenter(self):
        return self.center

    def getLength(self):
        return self.length

    def getHeight(self):
        return self.height

    def getSegmentLength(self):
        return self.length - (2.0*self.radius)

    def getRadius(self):
        return self.radius
