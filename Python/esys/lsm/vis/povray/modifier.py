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
Defines classes which relate to I{modifying} renderable objects.
"""

from esys.lsm.vis import core
from .color        import RgbColor

class Modifiable(core.Modifiable):
    def __init__(self):
        core.Modifiable.__init__(self)

    def apply(self, modifier):
        if (isinstance(modifier, RgbColor)):
            modifier = Pigment(modifier)
        core.Modifiable.apply(self, modifier)

    def applySingle(self, modifier):
        if (isinstance(modifier, RgbColor)):
            modifier = Pigment(modifier)
        core.Modifiable.applySingle(self, modifier)

    def writeSdl(self, f):
        for modifier in self.getModifierList():
            modifier.writeSdl(f)

class Modifier(core.Modifier):
    def __init__(self):
        pass

    def writeSdl(self, f):
        core.raiseNotImplementedError("Not implemented.")

class Pigment(Modifier):
    def __init__(self, pigment, transmit=0.0):
        self.pigment  = pigment
        self.transmit = transmit

    def writeSdl(self, f):
        f.write("pigment {")
        self.pigment.writeSdl(f)
        if (self.transmit > 0.0):
            f.write(" image_map { transmit ")
            f.write(str(self.transmit) + "}")
        f.write("}")

class Finish(Modifier):
    def __init__(
        self,
        ambient=0.2,
        diffuse=0.5,
        reflection=0.25,
        specular=0.4,
        roughness=0.01
    ):
        self.ambient    = ambient
        self.diffuse    = diffuse
        self.reflection = reflection
        self.specular   = specular
        self.roughness  = roughness

    def writeSdl(self, f):
        f.write("finish {")
        f.write(" ambient ")
        f.write(str(self.ambient))
        f.write(" diffuse ")
        f.write(str(self.diffuse))
        f.write(" reflection ")
        f.write(str(self.reflection))
        f.write(" specular ")
        f.write(str(self.specular))
        f.write(" roughness ")
        f.write(str(self.roughness))
        f.write("}")

class Orientation(Modifier):
    def __init__(self, matrix=None):
        self.matrix = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0], [0.0,0.0,0.0]]
        if (matrix != None):
            for i in range(0, matrix.getNumRows()):
                for j in range (0, matrix.getNumCols()):
                    self.matrix[j][i] = matrix(i,j)

    def writeSdl(self, f):
        f.write("matrix <")
        f.write("  {0:s},".format(",".join(map(str,self.matrix[0]))))
        f.write("  {0:s},".format(",".join(map(str,self.matrix[1]))))
        f.write("  {0:s},".format(",".join(map(str,self.matrix[2]))))
        f.write("  {0:s}".format(",".join(map(str,self.matrix[3]))))
        f.write(">")

class Scale(Modifier):
    def __init__(self, scale=1.0):
        if (hasattr(scale,"__iter__")):
            self.scale = tuple(scale)
        else:
            self.scale = (scale,scale,scale)

    def writeSdl(self, f):
        f.write(" scale <{0:s}> ".format(",".join(map(str,self.scale))))

class Pattern:
    def __init__(self):
        pass

class Checker(Pattern,Modifiable):
    def __init__(self, color0=None, color1=None):
        Pattern.__init__(self)
        Modifiable.__init__(self)
        self.color0 = color0
        self.color1 = color1

    def writeSdl(self, f):
        f.write("checker ")
        if (self.color0 != None):
            self.color0.writeSdl(f)
            f.write(" ")
        if (self.color1 != None):
            self.color1.writeSdl(f)
        Modifiable.writeSdl(self, f)



