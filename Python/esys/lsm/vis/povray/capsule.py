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
from esys.lsm.vis import core
from .modifier     import Modifiable
from .csg          import Union
from .box          import Box
from .cylinder     import Disk

class Capsule(core.Capsule, Modifiable):
    """
    Objects of this class represent capsules, a capsule capped with
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
        core.Capsule.__init__(self,radius,endPt1,endPt2,segEndPt1,segEndPt2)
        Modifiable.__init__(self)

    def writeSdl(self, f):
        radius = self.getRadius()
        f.write("\nsphere_sweep {\nlinear_spline\n2\n")
        f.write("<{0:s},0,0> {1:s}\n".format(-self.getSegmentLength()*0.5,radius))
        f.write("<{0:s},0,0> {1:s}\n".format(self.getSegmentLength()*0.5,radius))
        Modifiable.writeSdl(self, f)
        f.write("\ntranslate <{0:s},{1:s},{2:s}>\n".format(*tuple(self.getCenter())))
        f.write("}")

class CapsuleDisk(core.CapsuleDisk, Modifiable):
    """
    Objects of this class represent a cylinder with a capsule cross-section.
    """
    def __init__(
        self,
        center,
        length,
        radius,
        height = None
    ):
        core.CapsuleDisk.__init__(self,center,length,radius,height)
        Modifiable.__init__(self)

    def writeSdl(self, f):
        csg = Union()
        halfSegLength = 0.5*self.getSegmentLength()
        box = \
          Box(
            core.Vec3(-halfSegLength, -self.getRadius(), -0.5*self.getHeight()),
            core.Vec3(halfSegLength, self.getRadius(), 0.5*self.getHeight())
          )
        csg.append(box)
        csg.append(
          Disk(
            center = core.Vec3(-halfSegLength, 0, 0),
            radius = self.getRadius(),
            height = self.getHeight(),
            direction = (0,0,1)
          )
        )
        csg.append(
          Disk(
            center = core.Vec3(halfSegLength, 0, 0),
            radius = self.getRadius(),
            height = self.getHeight(),
            direction = (0,0,1)
          )
        )
        csg.writeBegin(f)
        csg.writeObjects(f)
        Modifiable.writeSdl(self,f)
        f.write("\ntranslate <{0:s},{1:s},{2:s}>\n".format(*tuple(self.getCenter())))
        csg.writeEnd(f)
