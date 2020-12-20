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

"""Defines the L{Cylinder} base class"""
from __future__   import division
from esys.lsm.vis import core
from .modifier     import Modifiable

class Cylinder(core.Cylinder, Modifiable):
    def __init__(self, endPt1, endPt2, radius):
        core.Cylinder.__init__(self, endPt1, endPt2, radius)
        Modifiable.__init__(self)

    def writeSdl(self, f):
        f.write("cylinder {")
        f.write(
            "\n<"
            +
            ",".join(
              map(
                str,
                core.Vec3(self.getEndPt1())-core.Vec3(self.getCenter())
              )
            )
            +
            "> "
        )
        f.write(
            "\n<"
            +
            ",".join(
              map(
                str,
                core.Vec3(self.getEndPt2())-core.Vec3(self.getCenter())
              )
            )
            +
            "> "
        )
        f.write(str(self.getRadius()))
        f.write(" ")
        Modifiable.writeSdl(self, f)
        f.write(" translate <")
        f.write(",".join(map(str,self.getCenter())))
        f.write(">")
        f.write("}")

class Disk(Cylinder):
    """
    Objects of this class represent cylinders.
    """
    def __init__(self, center, radius, height=None, direction=None):
        if (height == None):
           height = radius*0.01
        if (direction == None):
            direction = core.Vec3(0,0,1)
        direction = core.Vec3(direction)
        nrm = direction.norm()
        if (nrm == 0.0):
            direction = core.Vec3(0,0,1)
            nrm = 1.0

        center = core.Vec3(center)
        direction = direction/nrm
        Cylinder.__init__(
            self,
            endPt1 = center - (direction*(height*0.5)),
            endPt2 = center + (direction*(height*0.5)),
            radius = radius
        )
