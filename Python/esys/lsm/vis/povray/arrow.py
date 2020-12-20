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

"""Defines the L{Arrow} class"""
from __future__   import division
from esys.lsm.vis import core
from .cylinder     import Cylinder
from .cone         import Cone
from .csg          import Union
from .modifier     import Modifiable

class Arrow(core.Arrow, Modifiable):
    def __init__(self, tailPt, headPt, tailRadius, headRadius, headLength):
        core.Arrow.__init__(
            self,
            tailPt,
            headPt,
            tailRadius,
            headRadius,
            headLength
        )
        Modifiable.__init__(self)

    def writeSdl(self, f):
        tailPt = core.Vec3(self.getTailPt())
        tipPt = core.Vec3(self.getHeadPt())
        dir = tipPt - tailPt
        arrowLength = dir.norm()
        if (arrowLength > 0.0):
            u = Union()
            dir = dir/arrowLength
            basePt = tipPt - dir*self.getHeadLength()
            if (arrowLength > self.getHeadLength()):
                u.append(
                    Cylinder(
                        endPt1=tailPt,
                        endPt2=basePt,
                        radius=self.getTailRadius()
                    )
                )
            u.append(
                Cone(basePt=basePt, tipPt=tipPt, radius=self.getHeadRadius())
            )
            u.apply(self.getModifierList())
            u.writeSdl(f)
