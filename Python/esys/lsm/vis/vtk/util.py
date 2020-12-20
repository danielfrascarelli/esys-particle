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

"""
Common function and class definitions.
"""

from __future__   import division
import math
from esys.lsm.vis import core

def computeRotateWxyz(v1, v2):
    """
    Returns the angle (in degrees) and axis of rotation required
    to rotate v1 onto v2.
    @type v1: sequence of 3 floats
    @param v1: Vector to be rotated.
    @type v2: sequence of 3 floats
    @param v2: C{v1} rotated parallel to this vector.
    @rtype: tuple
    @return: C{(angle, axis)} rotation which rotates C{v1} parallel to C{v2}.
    """
    vv1 = core.Vec3(v1)
    vv2 = core.Vec3(v2)
    theta = (math.acos(vv2.dot(vv1)/(vv2.norm()*vv1.norm())))*(180.0/math.pi)
    rotAxis = vv2.cross(vv1)
    return (theta, rotAxis/(-(rotAxis.norm())))
