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

from esys.lsm          import LsmMpi, setVerbosity
from esys.lsm          import NRotBondedWallPrms
from esys.lsm.sim      import *
from esys.lsm.util     import *
from esys.lsm.geometry import *

"""
filter(..) in Python 3 constructs an iterator (as itertools.ifilter(..) in Python 2), but
a list in Python 2.
"""
if sys.version_info[0] > 2:
  iifilter = filter
else:
  from itertools import ifilter as iifilter

class MyExpSourcePrms(ExpSourcePrms):
    def getPosn(self, t):
        d = Vec3()
        for i in range(0, 3):
            d[i] = self.a[i]*math.exp(-(((t-self.t0[i])*self.b[i])**2))
        #print "Moving source by |" + str(d) + "|" + str(d.norm())
        return d

class MyCircularSourcePrms(SourcePrms):
    def __init__(self, posn, startTime=0.5, freq=0.05, radius=0.05):
        SourcePrms.__init__(self, posn)
        self.startTime = startTime
        self.freq      = freq
        self.radius    = radius

    def getPosn(self, t):
        d = Vec3(0,0,0)
        theta = self.freq*2.0*math.pi*(t-self.startTime)
        if ((t >= self.startTime) and (theta <= 2.0*math.pi)):
            d = \
                Vec3(
                    self.radius*math.cos(theta)-self.radius,
                    self.radius*math.sin(theta),
                    0
                )
            print("Moving source by |" + str(d) + "|=" + str(d.norm()))
        return d

class MyLinearSineSourcePrms(SourcePrms):
    def __init__(
        self,
        posn,
        startTime=0.5,
        freq=0.02,
        magnitude=[0, 0.10,0]
    ):
        SourcePrms.__init__(self, posn)
        self.startTime = startTime
        self.freq      = freq
        self.magnitude = magnitude

    def getPosn(self, t):
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

class PVisitor:
    """
    Objects of this class are used in conjunction with
    the waveProp.visitParticlesWithId method to collect model
    particle data.
    """
    def __init__(self):
        self.particleList = []

    def __iter__(self):
        return iter(self.particleList)

    def visitAParticle(self, particle):
        self.particleList.append(particle)

    def visitNRotSphere(self, p):
        self.visitAParticle(p)

    def visitRotSphere(self, p):
        self.visitAParticle(p)

  
def writeDisplacementData(
    idList,
    index,
    lsm,
    fileNamePrefix="displacement_"
):
    """
    Writes particle displacement data to file. Each line of the file is
    'px py pz dx dy dz' where p=(px,py,pz) is the particle position and
    d=(dx,dy,dz) is the current particle displacement (ie position relative
    to initial position).
    @param idList: list of particle id's.
    @param index: integer used to generate file name.
    @param lsm: a LSM object.
    @param fileNamePrefix: prefix of file where displacement data is saved.
    """
    #
    # Collect data for all particles with id's specified
    # in the idList.
    #
    visitor = PVisitor()
    lsm.visitParticlesWithId(idList, visitor)
    f = file(fileNamePrefix + "{0:d}.txt".format(index), "w")
    for p in visitor:
        f.write(
            str(p.getPosn()) +\
            " " + str(p.getPosn()-p.getInitialPosn()) + "\n"
        )
    f.close()
  
if (__name__ == "__main__"):
    # Create a 2d cubic close-packing of particles.
    radius = 1.0
    (nx,ny,nz) = (160, 160, 1)
    particles = CubicBlock([nx,ny,nz], radius)
    bBox = particles.getParticleBBox()
    centrePt = (bBox.getMinPt() + bBox.getMaxPt())*0.5
    
    #
    # Calculate bonds between particles. The DistConnections
    # object creates connections for a pair of particles which
    # are less than a specified distance appart.
    #
    connections = DistConnections(0.25, 0, particles)
    
    #
    # Generate lots of debug output by setting verbosity to True.
    #
    #setVerbosity(True)
    
    #
    # Create the wave propagation model,
    # Two worker processes, approx half the
    # the particles on one worker, half the particles
    # on the other worker, mpiDimList=[2,1,1] implies
    # splitting the domain in the 0 coordinate (x-coordinate).
    #
    waveProp = \
        WavePropagation(
            domainBox = bBox,
            do2d = (nz==1),
            numWorkerProcesses = 2,
            mpiDimList = [2,1,1],
            timeStepSize=0.04
        )
    
    #
    # Tag outer boundary particles so they can be bonded to
    # fixed walls.
    #
    tag = 1
    for p in \
      iifilter(
          lambda x:\
              abs(
                x.getPosn()[0]-x.getRadius()-bBox.getMinPt()[0]
              ) < 0.01,
          particles
      ):
      p.setTag(tag)
    tag += 1
    for p in \
      iifilter(
          lambda x:\
              abs(
                x.getPosn()[1]-x.getRadius()-bBox.getMinPt()[1]
              ) < 0.01,
          particles
      ):
      p.setTag(tag)
    tag += 1
    for p in \
      iifilter(
          lambda x:\
              abs(
                x.getPosn()[0]+x.getRadius()-bBox.getMaxPt()[0]
              ) < 0.01,
          particles
      ):
      p.setTag(tag)
    tag += 1
    for p in \
      iifilter(
          lambda x:\
              abs(
                x.getPosn()[1]+x.getRadius()-bBox.getMaxPt()[1]
              ) < 0.01,
          particles
      ):
      p.setTag(tag)
    
    #
    # Create the model particles
    #
    waveProp.createParticles(particles)
    
    #
    # Create the linear elastic bonds between particles.
    # The {0:1.0} argument is a dictionary with a
    # single (key=0, value=1.0) entry. All connections with
    # tag=0 will have a corresponding linear-elastic-bond created with
    # spring constant of 1.0.
    #
    waveProp.createBonds(connections, {0:1.0})
    
    #
    # Create the source disturbance which generates the wave.
    # A single particle is displaced over a small distance.
    # The source is created in the centre of the particle model.
    #
    approxSourcePosn = centrePt
    waveProp.createSources(MyLinearSineSourcePrms(approxSourcePosn))
    sourcePosn = waveProp.sourceList[0].getInitialPosn()
    print("Source posn = " + str(sourcePosn))
    
    #
    # Create the walls and the elastic bonds between walls
    # and the tagged particles.
    #
    wallBondSpringK = 1.0
    waveProp.createWall(
        NRotBondedWallPrms(
          wallBondSpringK, # spring constant
          bBox.getMinPt(), # plane/wall postition
          Vec3(1, 0, 0),   # plane/wall normal
          1                # particles with this tag get bonded to the wall
        )
    )
    waveProp.createWall(
        NRotBondedWallPrms(
          wallBondSpringK,
          bBox.getMinPt(),
          Vec3(0, 1, 0),
          2
        )
    )
    waveProp.createWall(
        NRotBondedWallPrms(
          wallBondSpringK,
          bBox.getMaxPt(),
          Vec3(-1, 0, 0),
          3
        )
    )
    waveProp.createWall(
        NRotBondedWallPrms(
          wallBondSpringK,
          bBox.getMaxPt(),
          Vec3(0, -1, 0),
          4
        )
    )
    
    #
    # Create a line of seismographs through the particle block
    #
    pt1 = Vec3(bBox.getMinPt())
    pt2 = Vec3(sourcePosn)
    
    numSeismos = 20
    diff = pt2-pt1
    interSeismoDistance = max([radius*2, diff.norm()/float(numSeismos)])
    incr = (diff/diff.norm())*interSeismoDistance
    seismographPosnList = []
    for i in range(0, numSeismos):
        seismographPosnList.append(pt1 + incr*float(i))
    waveProp.createSeismographGroup(
        seismographPosnList,
        "seismo_line_",
        sourcePosn
    )

    #
    # Run the model
    #
    numTimeSteps = 4000
    idList = [p.getId() for p in particles]
    j = 0
    t1 = None

    for i in range(0, numTimeSteps):
        if (t1 == None):
            t1 = time.time()
        waveProp.runTimeStep()
        if ((i % 5) == 0):
            t2 = time.time()
            waveProp.saveSeismoData()
            t3 = time.time()
            print("t = " + str(waveProp.getTime()) + ", step number = " + str(i) \
                + \
                ", seismo data output time = " \
                + \
                str(t3-t2) + " sec")
        if ((i % 50) == 0):
            t2 = time.time()
            writeDisplacementData(idList, j, waveProp)
            j += 1
            t3 = time.time()
            print("t = " + str(waveProp.getTime()) + ", step num = " + str(i) \
                + \
                ", run time = " + str(t2-t1) + " sec, displ data time = " \
                + \
                str(t3-t2) + " sec")
            t1 = None

    waveProp.writeReorderedRecordSectionData()
