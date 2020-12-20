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
Module for running a suite of gravity benchmark simulations.
Defines L{GravitySim} class whose instances can be used to run a
single simple gravitational simulation.
"""
from __future__          import division
from esys.lsm            import *
from esys.lsm.util       import *
from esys.lsm.geometry   import *
from esys.lsm.benchmarks import Util

from   time import time
import math
import sys

class GravitySim(Util.BenchSim):
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
        #
        # Set the time step size so that particles don't
        # move further than the verletDist.
        # This is so the sim doesn't need to update neighbour
        # lists.
        #
        self.gravAcceleration = Vec3(0,0,-9.8)
        self.logger.info("Initialising model")
        displDist = 2.0
        a = self.gravAcceleration.norm()
        timeStepSize = math.sqrt(displDist/a)/(2.0*numTimeSteps)
        self.initVerletModel(
          particleType = "NRotSphere",
          gridSpacing  = 4.0,
          verletDist = displDist
        )
        self.setTimeStepSize(dt = timeStepSize)

    def getSimTime(self):
        return self.getNumTimeSteps()*self.getTimeStepSize()

    def createParticles(self):
        """
        Creates a regularly packed block of particles within
        the model.
        """
        #
        # Create a block of regular packed particles.
        #
        blk = Util.getSimpleBlock(self.numParticles, radius=1.0)

        #
        # Now determine how far particles will displace under
        # gravitational acceleration and set the rectangular
        # simulation domain accordingly.
        #
        blkBox = blk.getParticleBBox()
        maxDisplacement = \
          self.gravAcceleration*self.getSimTime()*self.getSimTime()

        self.logger.info("Max displacement = {0:f}".format(maxDisplacement.norm()))

        absDisplacement = Vec3([abs(maxDisplacement[i]) for i in range(0,3)])
        domainBox = \
          BoundingBox(
            blkBox.getMinPt() - absDisplacement,
            blkBox.getMaxPt() + absDisplacement
          )
        self.logger.info("Setting spatial domain.")
        self.setSpatialDomain(domainBox)

        #
        # Finally create the model particles.
        #
        self.logger.info("Creating model particles.")
        LsmMpi.createParticles(self,blk)

    def createGravity(self):
        """
        Creates the gravitational body-force within the model.
        """
        self.logger.info("Creating gravity interaction group.")
        self.createInteractionGroup(
            GravityPrms(
              name="9.8m/s",
              acceleration = self.gravAcceleration
            )
        )

    def runSimulation(self):
        self.createParticles()
        self.createGravity()
        self.doTimedRun()
        
def runSimulations(
    numTimeSteps,
    numParticlesList,
    outputFileName,
    spawnCmdLineList,
    verboseLsm
):
    suite = Util.SimSuite(GravitySim)
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
