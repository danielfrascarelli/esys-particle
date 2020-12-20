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

"""Defines the L{Sphere} base class"""
from .exception import raiseNotImplemented

class Sphere(object):
    """
    Objects of this class represent spheres.
    """
    def __init__(self, center, radius):
        """
        Initialises sphere with center-coordinate and radius.
        """
        self.center = center
        self.radius = radius

    def getCenter(self):
        """
        Returns the coordinate of the center of this sphere.
        @return: Center coordinate of this sphere.
        """
        return self.center

    def getRadius(self):
        """
        Returns the radius of this sphere.
        @return: Radius of this sphere.
        """
        return self.radius
