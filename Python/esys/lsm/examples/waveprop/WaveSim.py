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

from __future__        import division, print_function
from sets              import Set
from itertools         import chain
import time
import sys

from esys.lsm          import LsmMpi, setVerbosity, NRotBondedWallPrms
from esys.lsm.sim      import *
from esys.lsm.util     import *
from esys.lsm.geometry import *

from .OptionParser     import OptionParser

"""
filter(..) in Python 3 constructs an iterator (as itertools.ifilter(..) in Python 2), but
a list in Python 2.
"""
if sys.version_info[0] > 2:
  iifilter = filter
else:
  from itertools import ifilter as iifilter

class MyLinearSineSourcePrms(SourcePrms):
    """
    Objects define linear spatial trajectory (sinusoidal
    time-frequency) for a particle-source-disturbance.
    """
    def __init__(self, posn, startTime=0.5, freq=0.02, magnitude=[0, 0.10, 0]):
        """
        Initialises trajectory parameters.
        @type posn: L{Vec3}
        @param posn: The approx location of the source-disturbance.
        @type startTime: float
        @param startTime: The time at which the disturbance begins.
        @type freq: float
        @param freq: The: frequency of the sinusoidal disturbance.
        @type magnitude: L{Vec3}
        @param magnitude: The direction and maximum displacement of
                         the source.
        """
        SourcePrms.__init__(self, posn)
        self.startTime = startTime
        self.freq      = freq
        self.magnitude = magnitude

    def getPosn(self, t):
        """
        Returns the position of the disturbance at the specified time.
        @type t: float
        @param t: Time for which position is returned.
        @rtype: L{Vec3}
        @return: Position of disturbance at time t.
        """
        d = Vec3(0,0,0)
        theta = self.freq*2.0*math.pi*(t-self.startTime)
        if ((t >= self.startTime) and (theta <= 1.0*math.pi)):
            d = \
                Vec3(
                    self.magnitude[0]*math.sin(theta),
                    self.magnitude[1]*math.sin(theta),
                    self.magnitude[2]*math.sin(theta)
                )
            print("Moving source by |" + str(d) + "|=" + str(d.norm()))
        return d

#
# Define particle tag values for particles which are to be
# elastically bonded to fixed walls
#
(LEFT_TAG, RIGHT_TAG, TOP_TAG, BOTTOM_TAG, FRONT_TAG, BACK_TAG) =\
    (101, 102, 103, 104, 105, 106)

UPPER_BOND_TAG = 10
LOWER_BOND_TAG = 11

class WaveSim:
    """
    Wrapper class for 2D and 3D wave propagation simulation.
    """
    def __init__(self, options=None):
        """
        Initialises simulation parameters.

        The C{options} argument, is an object with the following
        attributes:
            - C{radius} - Particle radius (float)
            - C{sourceDepth} - Relative depth of source disturbance (float in
              [0,1]).
            - C{sourceFrequency} - Frequency of sinusoidal source disturbance
              (float)
            - C{sourceMaxDisplacement} - Maximum relative displacement of source
              disturbance (list of 3 floats)
            - C{upperSpringK} - Spring constant for upper elastic medium
            - C{lowerSpringK} - Spring constant for lower elastic medium
            - C{upperMediumDepth} - Depth of upper elastic medium (float in
              [0,1])
            - C{particlesPerDim} - Size of particle block (list of 3 ints)
            - C{particleDataIncr} - Particle data is saved every
              C{particleDataIncr} time steps (int)
            - C{seismoDataIncr} - Seismograph data is saved every
              C{seismoDataIncr} time steps (int)
            - C{numWorkerProcesses} - Number of MPI worker processes (int)
            - C{mpiDimList} - Spatial division of domain amongst MPI workers
              (list of 3 ints)
            - C{timeStepSize} - Time step size for explicit integration (float)
            - C{maxNumTimeSteps} - Number of time steps to run the simulation
              (float)
            - C{verbosity} - If True output lots of LSM debugging info (bool)
        
        @type options: object
        @param options: An object with attributes as specified above.
        """
        if (options != None):
            self.radius          = options.radius
            self.sourceDepth     = options.sourceDepth
            self.sourceFrequency = options.sourceFrequency
            self.sourceMaxDisplacement=options.sourceMaxDisplacement
            self.upperSpringK    = options.upperSpringK
            self.lowerSpringK    = options.lowerSpringK
            self.upperMediumDepth= options.upperMediumDepth
            self.particlesPerDim = options.particlesPerDim
            self.particleDataIncr= options.particleDataIncr
            self.seismoDataIncr  = options.seismoDataIncr
            self.numWorkerProcesses = options.numWorkerProcesses
            self.mpiDimList      = options.mpiDimList
            self.timeStepSize    = options.timeStepSize
            self.maxNumTimeSteps = options.maxNumTimeSteps
            self.verbosity       = options.verbosity
        else:
            self.verbosity = False
            self.radius = 1.0
            self.sourceDepth = 0.5
            self.sourceFrequency=0.02
            self.sourceMaxDisplacement=[0.1,0.1,0.0]
            self.upperSpringK = 1.0
            self.lowerSpringK = 1.0
            self.upperMediumDepth = 0.25
            self.particlesPerDim = [160, 160, 1]
            self.particleDataIncr= 100
            self.seismoDataIncr  = 10
            self.numWorkerProcesses = 2
            self.mpiDimList = [2, 1, 1]
            self.timeStepSize = 0.05
            self.maxNumTimeSteps = 2000

        self.wallBondSpringK = 1.0
        self.particleBBox  = None
        self.particleBlock = None
        self.lsmWaveSim    = None

    def is3d(self):
        """
        Returns True if this is a 3D simulation (ie whether
        (self.particlesPerDim[2] > 1).
        @return: True if this is a 3D simulation, False if it is 2D.
        """
        return (self.particlesPerDim[2] > 1)

    def is2d(self):
        """
        Returns True if this is a 3D simulation (ie (not
        self.is3d()).
        @return: True if this is a 2D simulation, False if it is 3D.
        """
        return (not self.is3d())

    def createParticleBlock(self):
        """
        Creates a cubic close-packing of particles.
        @return: Collection of L{SimpleSphere} objects.
        """
        return CubicBlock(self.particlesPerDim, self.radius)

    def getParticleBlock(self):
        """
        Returns the collection/sequence of particles which represent an
        elastic block.
        @return: Collection of L{SimpleSphere} objects.
        """
        if (self.particleBlock == None):
            self.particleBlock = self.createParticleBlock()
        return self.particleBlock

    def getParticleBBox(self):
        """
        Returns bounding box of the elastic block of particles.
        @rtype: L{BoundingBox}
        @return: Axis aligned bounding box.
        """
        if (self.particleBBox == None):
            self.particleBBox = self.getParticleBlock().getParticleBBox()
        return self.particleBBox

    def tagBoundaryParticles(self):
        """
        Tags outer boundary particles so they can be bonded to
        fixed walls.
        """
        bBox = self.getParticleBBox()
        particles = self.getParticleBlock()

        distTol = 0.1*self.radius
        
        for p in \
          iifilter(
              lambda x:\
                  abs(
                    x.getPosn()[0]-x.getRadius()-bBox.getMinPt()[0]
                  ) < distTol,
              particles
          ):
          p.setTag(LEFT_TAG)

        for p in \
          iifilter(
              lambda x:\
                  abs(
                    x.getPosn()[1]-x.getRadius()-bBox.getMinPt()[1]
                  ) < distTol,
              particles
          ):
          p.setTag(BOTTOM_TAG)

        for p in \
          iifilter(
              lambda x:\
                  abs(
                    x.getPosn()[0]+x.getRadius()-bBox.getMaxPt()[0]
                  ) < distTol,
              particles
          ):
          p.setTag(RIGHT_TAG)

        for p in \
          iifilter(
              lambda x:\
                  abs(
                    x.getPosn()[1]+x.getRadius()-bBox.getMaxPt()[1]
                  ) < distTol,
              particles
          ):
          p.setTag(TOP_TAG)

        if (self.is3d()):
            for p in \
                iifilter(
                    lambda x:\
                        abs(
                          x.getPosn()[2]-x.getRadius()-bBox.getMinPt()[2]
                        ) < distTol,
                    particles
                ):
                p.setTag(BACK_TAG)

            for p in \
                iifilter(
                    lambda x:\
                        abs(
                          x.getPosn()[2]+x.getRadius()-bBox.getMaxPt()[2]
                        ) < distTol,
                    particles
                ):
                p.setTag(FRONT_TAG)

    def createLsmWaveSim(self):
        """
        Create the uninitialised wave propagation model ie the
        L{WavePropagation} object.
        """
        # Two worker processes imply that approx half the
        # the particles on one worker, half the particles
        # on the other worker, mpiDimList=[2,1,1] implies
        # splitting the domain in the 0 coordinate (x-coordinate).
        #
        self.lsmWaveSim = \
            WavePropagation(
                domainBox = self.getParticleBBox(),
                do2d = (self.is2d()),
                numWorkerProcesses = self.numWorkerProcesses,
                mpiDimList = self.mpiDimList,
                timeStepSize=self.timeStepSize
            )
        #
        # Generate lots of debug output by setting verbosity to True.
        #
        if (self.verbosity):
            setVerbosity(True)
        
    def getLsmWaveSim(self):
        if (self.lsmWaveSim == None):
            self.createLsmWaveSim()
        return self.lsmWaveSim

    def createParticles(self):
        """
        Create initial configuration of particles in the model
        """
        self.tagBoundaryParticles()
        self.getLsmWaveSim().createParticles(self.getParticleBlock())

    def createBonds(self):
        """
        Creates linear elastic bonds between particles. Two regions
        of bonds created with different elastic constants.
        """
        xzPlaneDepth =                                                 \
            self.getParticleBBox().getMaxPt()[1]                       \
            -                                                          \
            self.getParticleBBox().getSize()[1]*self.upperMediumDepth
        allParticles = Set(self.getParticleBlock())
        upperParticles = \
            Set(
                iifilter(
                    lambda p: p.getPosn()[1] >= xzPlaneDepth,
                    self.getParticleBlock()
                )
            )
        lowerParticles = allParticles.difference(upperParticles)

        #
        # The DistConnections object creates connections
        # for a pair of particles which are less than a
        # specified distance (self.radius/4.0) appart.
        #
        connections = DistConnections(self.radius/4.0)
        connections.addParticles(upperParticles, UPPER_BOND_TAG)
        connections.addParticles(lowerParticles, LOWER_BOND_TAG)

        #
        # Create the linear elastic bonds between particles.
        # The {TAG1:springK1, TAG2:springK2} argument is a dictionary
        # with a two (key, value) entries. All connections with
        # a tag of UPPER_BOND_TAG will have a corresponding
        # linear-elastic-bond created with spring constant of
        # self.upperSpringTag.
        #
        self.getLsmWaveSim().createBonds(
            connections,
            {
              UPPER_BOND_TAG:self.upperSpringK,
              LOWER_BOND_TAG:self.lowerSpringK
            }
        )

    def createSource(self):
        """
        Create the source disturbance which generates the wave.
        A single particle is displaced over a small distance.
        The source is created centred in the x and z coords,
        at a specified y coord depth.
        """
        bBox = self.getParticleBBox()
        size = bBox.getSize()
        centrePt = (bBox.getMinPt() + bBox.getMaxPt())/2.0
        approxSourcePosn = centrePt
        approxSourcePosn[1] = bBox.getMaxPt()[1] - self.sourceDepth*size[1]

        self.getLsmWaveSim().createSources(
            MyLinearSineSourcePrms(
              posn=approxSourcePosn,
              freq=self.sourceFrequency,
              magnitude=self.sourceMaxDisplacement
            )
        )
        sourcePosn = self.getLsmWaveSim().sourceList[0].getInitialPosn()
        print("Source posn = " + str(sourcePosn))

    def createBoundaryWalls(self):
        """
        Create the walls and the elastic bonds between walls
        and the tagged particles, leave the maximum y side as
        a free surface (no elastic wall).
        """
        bBox = self.getParticleBBox()
        waveProp = self.getLsmWaveSim()
        waveProp.createWall(
            NRotBondedWallPrms(
              self.wallBondSpringK, # spring constant
              bBox.getMinPt(),      # plane/wall postition
              Vec3(1, 0, 0),        # plane/wall normal
              LEFT_TAG          # particles with this tag get bonded to the wall
            )
        )
        waveProp.createWall(
            NRotBondedWallPrms(
              self.wallBondSpringK,
              bBox.getMinPt(),
              Vec3(0, 1, 0),
              BOTTOM_TAG
            )
        )
        waveProp.createWall(
            NRotBondedWallPrms(
              self.wallBondSpringK,
              bBox.getMaxPt(),
              Vec3(-1, 0, 0),
              RIGHT_TAG
            )
        )

        if (self.is3d()):
            waveProp.createWall(
                NRotBondedWallPrms(
                  self.wallBondSpringK,
                  bBox.getMinPt(),
                  Vec3(0, 0, 1),
                  BACK_TAG
                )
            )
            waveProp.createWall(
                NRotBondedWallPrms(
                  self.wallBondSpringK,
                  bBox.getMaxPt(),
                  Vec3(0, 0, -1),
                  FRONT_TAG
                )
            )

    def createLineOfSeismos(
        self,
        pt1,
        pt2,
        srcPt,
        numSeismos,
        fileNamePrefix
    ):
        """
        Create a line of seismographs through the particle block,
        between two specified points.
        @type pt1: L{Vec3}
        @param pt1: End point of line
        @type pt2: L{Vec3}
        @param pt2: other end-point of line
        @type srcPt: L{Vec3}
        @param srcPt: location of source point-disturbance.
        """

        diff = pt2-pt1
        #
        # Don't place seismos any closer together than two times
        # particle radius.
        #
        interSeismoDistance = \
            max(
                [
                  self.radius*2,
                  diff.norm()/float(numSeismos)
                ]
            )
        incr = (diff/diff.norm())*interSeismoDistance
        seismographPosnList = []
        for i in range(0, numSeismos):
            seismographPosnList.append(pt1 + incr*float(i))
        self.getLsmWaveSim().createSeismographGroup(
            seismographPosnList,
            fileNamePrefix,
            srcPt
        )

    def createSeismographs(self):
        """
        Creates lines of seismographs through the particle block
        """
        sourcePosn = self.getLsmWaveSim().sourceList[0].getInitialPosn()

        #
        # Create line of seismos from corner of block to the position
        # of the source disturbance.
        #
        bBox = self.getParticleBBox()
        pt1 = Vec3(bBox.getMaxPt())
        pt2 = Vec3(sourcePosn)
        self.createLineOfSeismos(pt1, pt2, sourcePosn, 20, "srcToCorner_")

        #
        # Create grid of numX by numZ seismos on the top surface (maximum y).
        #
        numX = 25
        numZ = 25
        if (self.is2d()):
            numZ = 1
        xDiff = bBox.getMaxPt()[0]-bBox.getMinPt()[0]
        xIncr = max(self.radius*2.0, xDiff/float(numX))
        zDiff = bBox.getMaxPt()[2]-bBox.getMinPt()[2]
        zIncr = max(self.radius*2.0, zDiff/float(numZ))

        xCoordList = [bBox.getMinPt()[0] + x*xIncr for x in range(0, numX)]
        zCoordList = [bBox.getMinPt()[2] + z*zIncr for z in range(0, numZ)]
        seismoPosnList = []
        y = bBox.getMaxPt()[1]
        for x in xCoordList:
            for z in zCoordList:
                seismoPosnList.append(Vec3(x, y, z))
        self.getLsmWaveSim().createSeismographGroup(
            seismoPosnList,
            "surfaceGrid_",
            sourcePosn
        )

    def getOutputParticleIdList(self):
        """
        Returns a list of particle ids, particles with id
        in this list have their position and displacement
        data written to file.
        @rtype: list of particle ids
        @return: List of particle ids.
        """
        return [p.getId() for p in self.getParticleBlock()]

    def runTimeSteps(self):
        """
        Runs the model for self.maxNumTimeSteps time steps.
        """
        numTimeSteps = self.maxNumTimeSteps
        idList = self.getOutputParticleIdList()
        waveProp = self.getLsmWaveSim()
        waveProp.setParticleDataIdList(idList)

        j = 0
        t1 = None

        for i in range(0, numTimeSteps):
            if (t1 == None):
                t1 = time.time()
            waveProp.runTimeStep()
            if ((i % self.seismoDataIncr) == 0):
                t2 = time.time()
                waveProp.saveSeismoData()
                t3 = time.time()
                print("t = " + str(waveProp.getTime()) +\
                    ", step number = " + str(i)  \
                    + \
                    ", seismo data output time = " \
                    + \
                    str(t3-t2) + " sec")
            if ((i % self.particleDataIncr) == 0):
                t2 = time.time()
                waveProp.writeParticleData(j, "particle_")
                j += 1
                t3 = time.time()
                print("t = " + str(waveProp.getTime()) + ", step num = " + str(i)\
                    + \
                    ", run time = " + str(t2-t1) + " sec, displ data time = " \
                    + \
                    str(t3-t2) + " sec")
                t1 = None
    
        waveProp.writeReorderedRecordSectionData()

    def runSim(self):
        """
        Sets up wave propagation model and executes simulation.
        """
        self.createLsmWaveSim()
        self.createParticles()
        self.createBoundaryWalls()
        self.createBonds()
        self.createSource()    
        self.createSeismographs()

        self.runTimeSteps()

def getOptionParser():
    """
    Returns the L{OptionParser} object useful for parsing command
    line options related to wave-propagation simulation.
    """
    return OptionParser()

__doc__=\
"""
Wave propagation simulation module. Contains WaveSim convenience class for
initialising LSM particle model and running wave propagation simulation.
This module may be run as a C{{__main__}} from the command line as
follows::
{0:s}
""".format("  " + str.replace(getOptionParser().format_help(), "\n", "\n  "))

if (__name__=="__main__"):
    parser = getOptionParser()
    (options, args) = parser.parse_args(sys.argv[1:])

    waveSim = WaveSim(options)
    waveSim.runSim()
