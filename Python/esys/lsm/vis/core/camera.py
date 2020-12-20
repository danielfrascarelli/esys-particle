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

"""Defines the L{Camera} base class"""
from .exception import raiseNotImplemented

from esys.lsm.util import Vec3

class Camera(object):
    """
    Objects of this class represent a camera.
    """
    def __init__(self):
        """
        Initialises camera object.
        """
        pass

    def getPosn(self):
        """
        Returns the position coordinate of this camera.
        @return: Position coordinate of this camera.
        """
        raiseNotImplemented()

    def setPosn(self, coord):
        """
        Sets the position coordinate of this camera.
        @type coord: iterable of 3 C{float} elements
        @param coord: New position of this camera.
        """
        raiseNotImplemented()

    def rotatePosn(self, axis, axisPt):
        """
        Rotates the camera position about the specified axis.
        @type axis: iterable of 3 C{float} elements
        @param axis: Axis of rotation and angle of rotation
        (C{ = axis.norm()} radians).
        @type axisPt: iterable of 3 C{float} elements
        @param axisPt: Axis of rotation is assumed to pass through this point.
        """
        self.setPosn(Vec3(self.getPosn()).rotate(Vec3(axis),Vec3(axisPt)))

    def rotate(self, axis, axisPt):
        """
        Rotates the camera position and look-at coordinate
        about the specified axis.
        @type axis: iterable of 3 C{float} elements
        @param axis: Axis of rotation and angle of rotation
        (C{ = axis.norm()} radians).
        @type axisPt: iterable of 3 C{float} elements
        @param axisPt: Axis of rotation is assumed to pass through this point.
        """
        self.setPosn(Vec3(self.getPosn()).rotate(Vec3(axis),Vec3(axisPt)))
        self.setLookAt(Vec3(self.getLookAt()).rotate(Vec3(axis),Vec3(axisPt)))

    def getLookAt(self):
        """
        Returns the coordinate at which this camera is looking.
        @return: The coordinate being looked-at by this camera.
        """
        raiseNotImplemented()

    def setLookAt(self, coord):
        """
        Sets the coordinate at which this camera is looking.
        @type coord: iterable of 3 C{float} elements
        @param coord: New coordinate to be looked-at by this camera.
        """
        raiseNotImplemented()

    def setZoom(self, factor):
        """
        Decrease the scale of the projection area by the specified factor.
        A value greater than 1 is a zoom-in, a value less than 1 is a zoom-out.
        @type factor: float
        @param factor: Zoom multiplier C{factor > 0.0}.
        """
        raiseNotImplemented()

    def getZoom(self):
        """
        Returns the zoom factor.
        @see: L{setZoom}
        @rtype: float
        @return: The zoom factor.
        """
        raiseNotImplemented()


