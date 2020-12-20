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
Defines the L{TriangulatedSurface} class.
"""
from esys.lsm.vis import core
from .modifier import Modifiable

class TriangulatedSurface(core.TriangulatedSurface, Modifiable):
    def __init__(self, nodeSequence, faceSequence):
        core.TriangulatedSurface.__init__(self, nodeSequence, faceSequence)
        Modifiable.__init__(self)

    def writeSdl(self, f):
        f.write("mesh2 {\n")
        f.write("  vertex_vectors {\n")
        f.write("    {0:s}".format(len(self.nodeSequence)))
        for node in self.nodeSequence:
            f.write(",\n    <{0:s},{1:s},{2:s}>".format(*tuple(node)))
        f.write("\n  }\n")
        f.write("  face_indices {\n")
        f.write("    {0:s}".format(len(self.faceSequence)))
        for face in self.faceSequence:
            f.write(",\n    <{0:s},{1:s},{2:s}>".format(*tuple(face)))
        f.write("\n  }\n")
        Modifiable.writeSdl(self, f)
        f.write("}")


