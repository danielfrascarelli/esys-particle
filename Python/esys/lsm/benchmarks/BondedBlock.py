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
Module for running a suite of bonded-block benchmark simulations.
Defines L{BondedBlockSim} class whose instances can be used to run a
single simple I{quivering} elastic block simulation.
"""
from __future__          import division
from esys.lsm            import *
from esys.lsm.util       import *
from esys.lsm.geometry   import *
from esys.lsm.benchmarks import Util

import random

class BondedBlockSim(Util.BenchSim):
    def __init__(
        self,
        numTimeSteps=10000,
        numParticles=1000,
        spawnExe = "",
        spawnArgList = None
    ):
        Util.BenchSim.__init__(
            self,
            numTimeSteps,
            numParticles,
            spawnExe,
            spawnArgList
        )
        self.normalK = 1.0
        self.radius = 1.0
        self.verletDist = self.radius/2.0
        self.logger.info("Initialising model")
        timeStepSize = 0.0001
        self.initVerletModel(
          particleType = "NRotSphere",
          gridSpacing  = 4.0,
          verletDist = self.verletDist
        )
        self.setTimeStepSize(dt=timeStepSize)

    def getRandomOffset(self):
        return Vec3([random.random() for i in range(0,3)])*(self.radius/20)

    def createParticlesAndBonds(self):
        """
        Creates a hexagonnaly close packed block of particles within
        the model along with linear elastic bonds connecting neighbouring
        particles.
        """
        #
        # Create a block of hcc particles.
        #
        blk = Util.getHexagBlock(self.numParticles, self.radius)

        #
        # Determine neighbouring particles
        #
        conns = DistConnections(0.01*self.radius, 0, blk)

        #
        # initialize RNG
        #
        random.seed(42)

        #
        # Perurb the initial positions of the particles.
        #
        for p in blk:
            p.setPosn(p.getPosn() + self.getRandomOffset())

        #
        # Set up the spatial domain.
        #
        blkBox = blk.getParticleBBox()
        domainBox = \
          BoundingBox(
            blkBox.getMinPt() - 2.0*self.radius,
            blkBox.getMaxPt() + 2.0*self.radius
          )
        self.logger.info("Setting spatial domain.")
        self.setSpatialDomain(domainBox)

        #
        # Create the model particles.
        #
        self.createParticles(blk)

        #
        # Create the linear elastic bonds between particles.
        #
        self.logger.info("Creating bonded interaction group...")
        self.createConnections(conns)
        self.createInteractionGroup(
            NRotBondPrms(
                tag=0,
                name="bonds",
                normalK=self.normalK,
                breakDistance=self.radius*10.0
            )
        )
        self.logger.info("Creating " + str(len(conns)) + " bonds...")

    def runSimulation(self):
        self.createParticlesAndBonds()
        self.doTimedRun()

def runSimulations(
    numTimeSteps,
    numParticlesList,
    outputFileName,
    spawnCmdLineList,
    verboseLsm
):
    suite = Util.SimSuite(BondedBlockSim)
    suite.createSims(
        numTimeSteps,
        numParticlesList,
        outputFileName,
        spawnCmdLineList,
        verboseLsm
    )
    suite.runSimulations()

if (__name__ == "__main__"):
    parser = Util.getOptionParser()
    (options, argList) = parser.parse_args()
    runSimulations(
      options.numTimeSteps,
      options.numParticlesList,
      options.outputFileName,
      options.spawnCmdLineList,
      options.verboseLsm
    )
