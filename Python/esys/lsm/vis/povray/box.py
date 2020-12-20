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

"""Defines the L{Box} and L{Cube} classes"""
from esys.lsm.vis import core
from .modifier     import Modifiable

class Box(core.Box, Modifiable):
    def __init__(self, minPt, maxPt):
        core.Box.__init__(self, minPt, maxPt)
        Modifiable.__init__(self)

    def writeSdl(self, f):
        halfSize = self.getSideLength()*0.5
        f.write("box {<" + ",".join(map(str, halfSize*-1.0)) + ">")
        f.write(",<" + ",".join(map(str, halfSize)) + "> ")
        Modifiable.writeSdl(self, f)
        f.write(" translate <")
        f.write(",".join(map(str,self.getCenter())))
        f.write(">")
        f.write("}")

class Cube(Box):
    """
    Box with all sides the same length.
    """
    def __init__(self, minPt, sideLength):
        """
        Constructs axis-aligned cube.
        @type minPt: sequence of 3 floats
        @param minPt: lower left back corner coordinate.
        @type sideLength: float
        @param sideLength: length of all sides.
        """
        Box.__init__(self, minPt, core.Vec3(minPt)+sideLength)
