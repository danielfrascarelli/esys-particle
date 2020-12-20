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
#!/bin/env mpipython
"""
Defines single class L{FixedLidWaveSim.FixedLidWaveSim} which extends
L{WaveSim.WaveSim}.
"""
from .WaveSim      import *
from .OptionParser import OptionParser

import time
import sys

class FixedLidWaveSim(WaveSim):
    """
    Extends L{WaveSim.WaveSim} class, overides the L{createBoundaryWalls}
    walls method simply by adding a "fixed lid" instead of having a
    free-surface.
    """
    def createBoundaryWalls(self):
        """
        Create the walls and the elastic bonds between walls
        and the tagged particles. Create an elastic wall at
        the surface (maximum y) in addition to walls on other
        sides of the particle block.
        """

        #
        # Create walls on all sides except the maximum y
        #
        WaveSim.createBoundaryWalls(self)

        #
        # Create the wall on the max y side
        #
        self.getLsmWaveSim().createWall(
            NRotBondedWallPrms(
              self.wallBondSpringK,
              self.getParticleBBox().getMaxPt(),
              Vec3(0, -1, 0),
              TOP_TAG
            )
        )

if (__name__=="__main__"):
    parser = OptionParser()
    (options, args) = parser.parse_args(sys.argv[1:])

    waveSim = FixedLidWaveSim(options)
    waveSim.runSim()
