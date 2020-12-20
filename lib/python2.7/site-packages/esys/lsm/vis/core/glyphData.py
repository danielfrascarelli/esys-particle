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
Defines L{GlyphData} base class.
"""
class GlyphData(object):
    """
    An object which is rendered as multiple glyphs (of the same type)
    at a collection of points in space.
    """
    def __init__(self, data, glyphExtractor):
        """
        Constructs the object.
        @type data: iterable
        @param data: The data which is to be rendered as multiple glyphs.
        @type glyphExtractor: 
        """
        self.data           = data
        self.glyphExtractor = glyphExtractor

    def getData(self):
        """
        Returns the C{data} object associated with this C{GlyphData}
        object.
        @rtype: iterable
        @return: The C{data} object which is used to generate glyphs.
        """
        return self.data

    def getGlyphExtractor(self):
        """
        Returns the object which is used to extract individual
        glyph data from each element of the C{self.getData()} iterable.
        """
        return self.glyphExtractor
