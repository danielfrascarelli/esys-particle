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

from __future__         import division
from .color             import Colors
from esys.lsm.vis.core import Vec3

def writeVec(f, vec):
    f.write("<")
    f.write(str.join(",", (list(map(str,vec)))))
    f.write(">")

class Background(object):
    def __init__(self, color):
        self.color = color

    def writeSdl(self, f):
        f.write("background {")
        self.color.writeSdl(f)
        f.write("}")

class Plane(object):
    def __init__(self, posn, normal):
        self.posn   = Vec3(posn)
        self.normal = Vec3(normal)
        self.normal = self.normal/(self.normal.norm())
        self.pigment = None

    def setPigment(self, pigment):
        self.pigment = pigment

    def writeSdl(self, f):
        d = self.posn.dot(self.normal)
        f.write("\nplane {")
        writeVec(f, self.normal)
        f.write(",{0:s} ".format(d))
        if (self.pigment != None):
            self.pigment.writeSdl(f)
        f.write("}")

class LightSource(object):
    def __init__(self, posn, color=Colors.White, type="shadowless"):
        self.posn  = posn
        self.color = color
        self.type  = type

    def getPosn(self):
        return self.posn

    def setPosn(self, posn):
        self.posn = posn

    def getColor(self):
        return self.color

    def setColor(self, color):
        self.color = color

    def writeSdl(self, f):
        f.write("light_source {\n<")
        f.write(str.join(",", (list(map(str,self.getPosn())))))
        f.write(">\n")
        self.getColor().writeSdl(f)
        f.write("\n")
        f.write(str(self.type))
        f.write("\n}")
