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

class GlyphData(core.GlyphData):
    def __init__(self, data, glyphExtractor):
        core.GlyphData.__init__(self, data, glyphExtractor)

    def addActor(self, addTo):
        extractor = self.getGlyphExtractor()
        for record in self.getData():
            extractor.addActor(record, addTo)
