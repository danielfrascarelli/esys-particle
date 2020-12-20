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

from esys.lsm.vis import core
from .modifier     import Finish, Pigment
from .color        import Colors

class GlyphData(core.GlyphData):
    def __init__(self, data, glyphExtractor):
        core.GlyphData.__init__(self, data, glyphExtractor)
        self.finish = Finish()
        self.pigment = Pigment(Colors.Red)

    def getFinish(self):
        return self.finish

    def setFinish(self, finish):
        self.finish = finish

    def getPigment(self):
        return self.pigment

    def setPigment(self, pigment):
        self.pigment = pigment

    def writeSdl(self, f):
        f.write("union {\n")
        for record in self.getData():
            self.getGlyphExtractor().writeSdl(f, record)
            f.write("\n")
        self.getFinish().writeSdl(f)
        f.write("\n")
        self.getPigment().writeSdl(f)
        f.write("\n}")
