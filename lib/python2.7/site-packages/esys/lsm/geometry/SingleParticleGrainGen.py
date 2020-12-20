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
Defines the L{SingleParticleGrainGen} class which extends
L{RndGrainGen}.
"""
from esys.lsm.util                import *
from esys.lsm.geometry.GeometryPy import *

class SingleParticleGrainGen(RndGrainGen):
    """
    Grain generator which can be used in conjunction with
    the L{GougeConfigPrms} and L{GrainRandomBoxPacker} classes.
    Generates grains which consist of a single particle.
    @see: getGrain
    """
    def __init__(self, minGRadius, maxGRadius):
        """
        Initialise generator.
        @type minGRadius: float
        @param minGRadius: minimum grain radius
        @type maxGRadius: float
        @param maxGRadius: maximum grain radius
        """
        RndGrainGen.__init__(
            self,
            minGRadius,
            maxGRadius,
            minGRadius,
            maxGRadius
        )

    def getGrain(self, p):
        """
        Returns a grain consisting of a single particle
        which is I{copy-constructed} from C{p}.
        @type p: L{SimpleSphere}
        @param p: a particle
        @rtype: L{Grain}
        @return: A grain consisting of a single particle.
        """
        g = Grain()
        g.createParticle(p)
        return g

