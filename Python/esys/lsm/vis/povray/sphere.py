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
from esys.lsm.vis import core
from .modifier     import Modifiable

class Sphere(core.Sphere, Modifiable):
    def __init__(self, center, radius):
        core.Sphere.__init__(self, center, radius)
        Modifiable.__init__(self)

    def writeSdl(self, f):
        f.write("sphere {<0,0,0> ")
        f.write(str(self.getRadius()))
        f.write("\n")
        Modifiable.writeSdl(self, f)
        f.write("\ntranslate <")
        f.write(",".join(map(str,self.getCenter())))
        f.write(">")
        f.write("}")
