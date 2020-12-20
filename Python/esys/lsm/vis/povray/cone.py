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

"""Defines the L{Cone} base class"""
from esys.lsm.vis import core
from .modifier     import Modifiable

class TruncatedCone(core.TruncatedCone, Modifiable):
    def __init__(self, endPt1, endPt2, radius1, radius2 = 0):
        core.TruncatedCone.__init__(
            self,
            endPt1,
            endPt2,
            radius1,
            radius2
        )
        Modifiable.__init__(self)

    def writeSdl(self, f):
        f.write("cone {")
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
        f.write(" {0:s}".format(self.getRadius1()))
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
        f.write(str(self.getRadius2()))
        f.write(" ")
        Modifiable.writeSdl(self, f)
        f.write("\ntranslate <")
        f.write(",".join(map(str,self.getCenter())))
        f.write(">\n")
        f.write("}")

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



