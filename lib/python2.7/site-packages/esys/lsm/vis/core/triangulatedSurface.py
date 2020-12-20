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
class TriangulatedSurface(object):
    def __init__(self, nodeSequence, faceSequence):
        """
        Initialise the mesh data.
        @type nodeSequence: iterable
        @param nodeSequence: The list of triangle vertices.
        @type faceSequence: iterable
        @param faceSequence: A list three-tuples, with each three-tuple
        representing a triangular facet. Each three-tuple element is
        an index into the nodeSequence array indicating a vertex of
        the triangular facet.
        """
        self.nodeSequence = list(nodeSequence)
        self.faceSequence = list(faceSequence)


