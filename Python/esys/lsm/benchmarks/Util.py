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
Common functionality between benchmarks.
"""

from __future__        import division, print_function
from time              import time
import math
import sys

from esys.lsm          import *
from esys.lsm.util     import *
from esys.lsm.geometry import *

"""
range(..) in Python 3 returns an iterable (as xrange() in Python 2), but
a list in Python 2.
"""
if sys.version_info[0] > 2:
  iRange = range
else:
  iRange = xrange

logger = Logging.getLogger("benchmarks.Util")

def getDimList(numParticles):
    ppd = int(math.ceil(math.pow(numParticles, (1.0/3.0))))
    return [ppd,ppd,ppd]

def getParticleCollection(numParticles, iterable):
    pIter = iter(iterable)
    pColl = ParticleCollection()
    for i in iRange(0, numParticles):
        pColl.createParticle(next(pIter))
    return pColl
    
def getSimpleBlock(numParticles, radius):
    return \
        getParticleCollection(
            numParticles,
            SimpleBlock(getDimList(numParticles), radius)
        )

def getHexagBlock(numParticles, radius):
    return \
        getParticleCollection(
            numParticles,
            HexagBlock(getDimList(numParticles), radius)
        )

class BenchSim(LsmMpi):
    def __init__(
        self,
        numTimeSteps=10000,
        numParticles=1000,
        spawnExe = "",
        spawnArgList = None
    ):
        self.logger = Logging.getLogger(str(self.__class__.__name__))
        if (spawnArgList == None):
            spawnArgList = []
        self.logger.info("Constructing LsmMpi")
        LsmMpi.__init__(self,1,[1,0,0], spawnExe, spawnArgList)
        self.logger.info("Setting number of time steps to {0:d}".format(numTimeSteps))
        self.setNumTimeSteps(numTimeSteps)
        self.logger.info("Initialising attributes")
        self.numParticles = numParticles
        self.startTime = 0
        self.stopTime  = 0
        self.runTime   = 0

    def doTimedRun(self):
        logger.info("Getting number of particles.")
        numParticles = self.getNumParticles()
        logger.info("Running simulation, {0:d} particles.".format(numParticles))
        self.startTime = time()
        self.run()
        self.stopTime = time()
        self.runTime = self.stopTime - self.startTime
        if (numParticles != self.getNumParticles()):
            raise \
                Exception(
                    "Mismatch between number of particles at start of" +\
                    " simulation ({0:d}) and number at end of simulation ({1:d})."\
                    .format(numParticles, self.getNumParticles())
                )

    def getRunTime(self):
        return self.runTime

class SimSuite:
    def __init__(self, benchSimClass):
        self.simClass = benchSimClass
        self.simList  = []
        self.outputFileName = None
        self.logger = Logging.getLogger("benchmarks.SimSuite")

    def __iter__(self):
        return iter(self.simList)

    def runSimulations(self):
        self.logger.info("Running simulations...")
        for sim in self.simList:
            sim.runSimulation()
            print("Run time for ",sim.getNumParticles(), " = ",sim.getRunTime())

        if (self.outputFileName != None):
            f = file(self.outputFileName, "w")
        else:
            f = sys.stdout
        for sim in self.simList:
            f.write(str(sim.getNumParticles()) + " ")
            f.write(str(sim.getRunTime()) + "\n")

    def run(self):
        self.runSimulations()


    def createSim(
      self,
      numTimeSteps,
      numParticles,
      spawnCmd,
      spawnCmdArgList
    ):
        return \
            self.simClass(
                numTimeSteps,
                numParticles,
                spawnCmd,
                spawnCmdArgList
            )

    def createSims(
        self,
        numTimeSteps,
        numParticlesList,
        outputFileName,
        spawnCmdLineList = None,
        verboseLsm = False
    ):
        setVerbosity(verboseLsm)
        self.outputFileName = outputFileName
        if ((spawnCmdLineList == None) or (len(spawnCmdLineList)==0)):
            spawnCmdLineList = [""]
        self.logger.info("Creating " + str(self.simClass) + " simulations...")
        self.simList = \
          [
            self.createSim(numTimeSteps,numParticles,spawnCmdLineList[0],spawnCmdLineList[1:])
            for numParticles in numParticlesList
          ]

def getOptionParser():
    parser = \
        OptParse.LogOptionParser(
            usage =\
              "usage: %prog [options]\n\n" +\
              "Runs a suite of benchmark simulations for varying number of particles. "
        )
    parser.add_option(
      "-n", "--num-time-steps",
      dest="numTimeSteps",
      type="int",
      metavar="N",
      default=10000,
      help=\
          "The number of time steps for which a simulation is run."+\
          " (default N=10000)"
    )
    parser.add_option(
      "-p", "--num-particles-list",
      dest="numParticlesList",
      type="int_list",
      metavar="L",
      default=[1000],
      help=\
          "A list specifying the number of particles in each simulation"+\
          " (default L=[1000]). One simulation is run for each element in"+\
          " this list."
    )
    parser.add_option(
      "-o", "--output-file-name",
      dest="outputFileName",
      type="string",
      metavar="F",
      default=None,
      help=\
          "Timing data for each simulation is written to this file"+\
          " (default F=None)."
    )
    parser.add_option(
      "-s", "--spawn-cmd-line-list",
      dest="spawnCmdLineList",
      type="string_list",
      metavar="C",
      default=None,
      help=\
          "A list specifying the command and arguments used to spawn MPI"+\
          " worker processes (default L=[1000])."
    )
    parser.add_option(
      "-v", "--verbose-lsm",
      dest="verboseLsm",
      action="store_true",
      default=False,
      help=\
          "Generates (lots of) debug output during simulation run" +\
          " (default no debug output)"
    )

    return parser
