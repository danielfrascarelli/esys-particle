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

from __future__    import division, print_function
from sets          import Set
import math
import weakref
import sys


from esys.lsm.util import *
from esys.lsm      import NRotBondPrms, DampingPrms, Runnable
import WavePropagationPy

"""
filter(..) in Python 3 constructs an iterator (as itertools.ifilter(..) in Python 2), but
a list in Python 2.
"""
if sys.version_info[0] > 2:
  iifilter = filter
else:
  from itertools import ifilter as iifilter

class SourcePrms:
    """
    Base class for point-source distubance which generates elastic wave.
    """
    def __init__(self, posn):
        """
        @param posn: The desired location for the source disturbance.
        """
        self.initialPosn = Vec3(posn)

    def getInitialPosn(self):
        """
        Returns the desired location of the source disturbance.
        """
        return self.initialPosn

    def getPosn(self, t):
        """
        Returns relative position of source disturbance for time t.
        @param t: time.
        """
        raise NotImplementedException()

class ExpSourcePrms(SourcePrms):
    """
    Describes (time, position) trajectory information for a source particle.
    """
    def __init__(
        self,
        posn = Vec3(0.0, 0.0, 0.0),
        a  = (0.20, 0.25, 0.00),
        b  = (0.50, 0.50, 0.50),
        t0 = (5.00, 4.00, 6.00)
    ):
        """
        Constructs and exponential (Gaussian in time) source.
        @type posn: L{Vec3}
        @param posn: coordinate location of the source.
        @type a: L{list}/L{tuple}
        @param a: sequence of 3 elements defining the magnitude/amplitude of
                 the (x,y,z) components.
        @type b: L{list}/{tuple}
        @param b: sequence of 3 elements defining the period of the (x,y,z)
                 components.
        @type t0: float
        @param t0: sequence of 3 elements defining the time at which the maximum
                 amplitude is achieved for each (x,y,z) component.
        """
        SourcePrms.__init__(self, posn)
        self.a           = Vec3(a)
        self.b           = Vec3(b)
        self.t0          = Vec3(t0)

    def getPosn(self, t):
        d = Vec3()
        for i in range(0, 3):
            d[i] = self.a[i]*math.exp(-((t-self.t0[i])/self.b[i])**2)
        return d

class CircularSourcePrms(SourcePrms):
    """
    Defines circular trajectory (time, position) for source disturbance.
    """
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

class WaveSource(Runnable):
    """
    Helper class which moves a source particle during each time step.
    """
    def __init__(self, prms, lsm):
        """
        Constructs the source Runnable.
        @param prms: object which governs the source particle trajectory.
                     This object is expected to define a getPosn(t) method
                     which returns a relative position for time t, and
                     also define a getInitialPosn() method which returns
                     the location of the source disturbance.
        @type lsm: L{esys.lsm.sim.wavePropagation.WavePropagation}
        @param lsm: a WavePropagation lattice solid model object.
                    A source particle in this model is moved according to
                    the trajectory specifed by the prms argument.
        """
        Runnable.__init__(self)
        self.prms = prms
        self.lsmProxy = weakref.proxy(lsm)
        self.particleId = \
            self.getLsm().findClosestParticle(
                self.getPrms().getInitialPosn()
            )
        self.initialPosn = \
            self.getLsm().getParticlePosn(self.getParticleId())

    def getLsm(self):
        return self.lsmProxy

    def getPrms(self):
        return self.prms

    def getParticleId(self):
        return self.particleId

    def getInitialPosn(self):
        return self.initialPosn

    def run(self):
        posn =                    \
            self.getInitialPosn() \
            +                     \
            self.getPrms().getPosn(self.getLsm().getTime())
        self.getLsm().moveParticleTo(self.getParticleId(), posn)

class Seismograph:
    """
    Helper class whose objects represent seismographs.
    """
    def __init__(self, posn, lsm):
        """
        @type posn: Three float squence
        @param posn: the approximate location of the seismograph.
        @param lsm: a lattice solid model object.
        """
        self.lsmProxy = weakref.proxy(lsm)
        self.particleId = \
            self.getLsm().findClosestParticle(posn)
        self.posn = self.getLsm().getParticlePosn(self.particleId)

    def getLsm(self):
        """
        Returns LSM object which is used to obtain displacement,
        velocity and acceleration data.
        """
        return self.lsmProxy

    def getParticleId(self):
        """
        Returns the id of the particle which is used to obtain
        seismograph data.
        """
        return self.particleId

    def getInitialPosn(self):
        """
        Returns the initial position of the of this seismograph,
        this position is just the particle initial position.
        """
        return self.posn

def cmpRecordSetElemList(elemList1, elemList2):
    """
    Compare-function used with the list.sort method to order
    record section data by distance-to-source.
    """
    c = cmp(elemList1[1],elemList2[1])
    if (c == 0):
        return cmp(elemList1[0],elemList2[0])
    return c

class SeismographGroup:
    """
    Objects of this class represent a collection of seismographs.
    """
    def __init__(self, seismoList, name, fileNamePrefix, sourcePosn):
        """
        @type seismoList: list
        @param seismoList: List of Seismograph objects.
        @type name: string
        @param name: Name of this group.
        @type fileNamePrefix: string
        @param fileNamePrefix: For all Seismograph objects, data is
               saved to files with this prefix.
        @type sourcePosn: L{esys.lsm.util.FoundationPy.Vec3}
        @param sourcePosn: Location of source disturbance, used to
                           calculate distance from source to seismograph.
        """
        self.name           = name
        self.fileNamePrefix = fileNamePrefix
        self.seismoList     = seismoList
        self.seismoIdSet    = Set([s.getParticleId() for s in self.seismoList])
        self.sourcePosn     = Vec3(sourcePosn)

    def inGroup(self, seismoData):
        """
        Returns whether the specfied SeismographData object is
        from a Seismograph in this group.
        @type seismoData: L{SeismographData}
        @param seismoData: Determine whether this data belongs to a
                           Seismograph object in this group.
        @return: True if seismoData came from a seismograph in this group.
        """
        return (seismoData.particleId in self.seismoIdSet)

    def getRecordSectionFileName(self):
        """
        Returns the name of the file in which record section data
        is saved.
        """
        return self.fileNamePrefix + "record_section.txt"
    
    def getSeismoDataFileName(self, posn):
        """
        Returns the name of a seismograph data file for the specified posn.
        @param posn: A Vec3 object specifying the position of the seismograph.
        @rtype string
        @return: File name for seismograph data at spatial coordinate posn.
        """
        fileName = self.fileNamePrefix + "{0:3.3f}_{1:3.3f}_{2:3.3f}".format(*(posn.toTuple()))
        fileName += ".txt"
        return fileName

    def initialiseSeismoDataFiles(self):
        """
        Creates/overwrites empty seismograph data files.
        """
        fileName = self.getRecordSectionFileName()
        f = file(fileName, "w")
        f.close()
        for seismo in self.seismoList:
            fileName = self.getSeismoDataFileName(seismo.getInitialPosn())
            f = file(fileName, "w")
            f.write(
              "# " + "{0:3.3f} {1:3.3f} {2:3.3f}"\
              .format(*(seismo.getInitialPosn().toTuple()))
            )
            f.close()

    def save(self, time, seismoDataList):
        """
        Saves seismograph data to files. Appends record section info
        as well as individual seismo data file.
        """
        for seismoData in iifilter(self.inGroup, seismoDataList):
            self.saveSeismoData(time, seismoData)
            self.saveRecordSectionData(time, seismoData)
    
    def saveSeismoData(self, time, seismoData):
        """
        Appends the specfied seismograph data to file. Each seismograph data
        record is appended to a different file. Each line in the file is
        of the form 't dx dy dz vx vy vz ax ay az' where d=(dx,dy,dz) is
        the displacement vector, v=(vx,vy,vz) is the velocity vector and
        a=(ax,ay,az) is the acceleration vector.
        @type time: float
        @param time: the time the data was recorded
        @type seismoData: L{SeismographData}
        @param seismoData: Data in this object is saved to file. The file
                           name is determined by the data's position.
        """
        fileName = \
            self.getSeismoDataFileName(
                seismoData.posn
            )
        f = file(fileName, "a")
        f.write(str(time))
        f.write(" " + str(seismoData.displacement))
        f.write(" " + str(seismoData.velocity))
        f.write(" " + str(seismoData.acceleration))
        f.write("\n")
        f.close()

    def saveRecordSectionData(
        self,
        time,
        seismoData
    ):
        """
        Appends record-set style seismograph data to file. Each seismograph data
        record is appended to a record-set file. Each line in the file is
        of the form 't d ux uy uz vx vy vz ax ay az' where t is the time, d is
        the distance from the seismograph to the sourcePosn d=(dx,dy,dz)
        is the displacement vector, v=(vx,vy,vz) is the velocity vector and
        a=(ax,ay,az) is the acceleration vector.
        @type time: float
        @param time: the time the data was recorded
        @type seismoData: L{SeismographData}
        @param seismoData: This data is saved to record the section file.
        """
        fileName = self.getRecordSectionFileName()
        f = file(fileName, "a")
        f.write(str(time))
        f.write(" " + str((seismoData.posn-self.sourcePosn).norm()))
        f.write(" " + str(seismoData.displacement))
        f.write(" " + str(seismoData.velocity))
        f.write(" " + str(seismoData.acceleration))
        f.write("\n")
        f.close()

    def writeReorderedRecordSectionData(self):
        """
        Re-orders record set data for more convenient plotting in gnuplot.
        Reads time-ordered record section data and re-orders according to
        distance-to-source.
        """
        f = file(self.getRecordSectionFileName(), "r")
        recordList = []
        for line in f.readlines():
            line = str.strip(line)
            if (len(line) > 0):
                recordList.append(list(map(float, str.split(line, " "))))
        f.close()
        #
        # Separate each seismo's data into a list
        #
        distDict = dict()
        for record in recordList:
            if (record[1] not in distDict):
                distDict[record[1]] = []
            distDict[record[1]].append(record)
        f = file(self.fileNamePrefix + "record_section_reordered.txt", "w")
        distKeyList = list(distDict.keys())
        distKeyList.sort()
        for dist in distKeyList:
            distDict[dist].sort(cmpRecordSetElemList)
            for recordElemList in distDict[dist]:
                f.write(str.join(" ", (list(map(str, recordElemList)))))
                f.write("\n")
            f.write("\n")
        f.close()

class SeismographGroupCollection:
    """
    Represents a collection/list of SeismographGroup objects. 
    """
    def __init__(self, lsm):
        """
        @type lsm: L{esys.lsm.sim.wavePropagation.WavePropagation}
        @param lsm: Seismo data is collected from this object.
        """
        self.lsmRef = weakref.ref(lsm)
        self.idSeismoDict = dict()
        self.seismoGroupList = []

    def getLsm(self):
        """
        Returns the LSM object associated with this collection.
        @rtype: L{esys.lsm.sim.wavePropagation.WavePropagation}
        """
        return self.lsmRef()

    def createSeismograph(self, posn):
        """
        Creates a seismograph at the specified position.
        @type posn: L{Vec3}
        @param posn: Approximate location of seismograph.
        @rtype: L{Seismograph}
        @return: Newly created Seismograph object.
        """
        return Seismograph(posn, self.getLsm())

    def getSeismograph(self, posn):
        """
        Returns seismograph for the given position. Only allows a
        unique instance of a seismograph for a particular spatial
        region.
        @rtype: L{Seismograph}
        """
        seismo = self.createSeismograph(posn)
        if (seismo.getParticleId() not in self.idSeismoDict):
            self.idSeismoDict[seismo.particleId] = seismo
        return self.idSeismoDict[seismo.particleId]

    def createGroup(self, posnIterable, fileNamePrefix, sourcePosn):
        """
        Creates multiple seismographs from a sequence of locations.
        @type posnIterable: sequence of L{Vec3}
        @param posnIterable: sequence of seismograph approximate locations.
        @type fileNamePrefix: string
        @param fileNamePrefix: prefix of file where seimo data is saved.
        @type sourcePosn: L{Vec3}
        @param sourcePosn: Location of source disturbance.
        """
        seismoList = [self.getSeismograph(p) for p in posnIterable]
        self.seismoGroupList.append(
            SeismographGroup(
                seismoList,
                str(len(self.seismoGroupList)),
                fileNamePrefix,
                sourcePosn
            )
        )
        self.seismoGroupList[-1].initialiseSeismoDataFiles()

    def getSeismographData(self):
        """
        Returns a sequence of SeismographData objects representing
        seismograph data for the current time step.
        """
        visitor = ParticleDataVisitor()
        self.getLsm().visitParticlesWithId(list(self.idSeismoDict.keys()), visitor)
        return visitor.seismoDataList

    def saveSeismoData(self):
        """
        Saves seismograph data for each group of seismographs.
        """

        #
        # Start by getting the data for all seismographs so
        # there is only one MPI gather call to the worker processes.
        #
        time = self.getLsm().getTime()
        seismoDataList = self.getSeismographData()

        #
        # Now each group writes the data to the appropriate files.
        #
        for seismoGroup in self.seismoGroupList:
            seismoGroup.save(time, seismoDataList)

    def writeReorderedRecordSectionData(self):
        """
        Writes reordered seismo data to file, data is ordered
        by distance-to-source then by time.
        """
        for seismoGroup in self.seismoGroupList:
            seismoGroup.writeReorderedRecordSectionData()

class SeismographData:
    """
    Objects of this class represent seismograph data at a particular
    instant in time.
    """
    def __init__(self, posn, displacement, velocity, acceleration, particleId):
        self.posn         = Vec3(posn)
        self.displacement = Vec3(displacement)
        self.velocity     = Vec3(velocity)
        self.acceleration = Vec3(acceleration)
        self.particleId   = particleId

class ParticleDataVisitor:
    """
    Helper class whose objects are used to collect
    seismograph data, from multiple locations,
    at a particular instant in time.
    Used in conjunction with the
    esys.lsm.sim.WavePropagation.WavePropagation.visitParticlesWithId
    method.
    """
    def __init__(self):
        self.seismoDataList = []

    def visitAParticle(self, p):
        """
        Converts the Particle data object p into a
        SeismographData object.
        @type p: L{esys.lsm.LsmPy.NRotSphere} or L{esys.lsm.LsmPy.RotSphere}
        @param p: Particle being visited.
        """
        self.seismoDataList.append(
            SeismographData(
                posn         = p.getInitialPosn(),
                displacement = p.getPosn() - p.getInitialPosn(),
                velocity     = p.getVelocity(),
                acceleration = p.getAcceleration(),
                particleId   = p.getId()
            )
        )

    def visitNRotSphere(self, p):
        self.visitAParticle(p)

    def visitRotSphere(self, p):
        self.visitAParticle(p)

class PVisitor:
    """
    Objects of this class are used in conjunction with
    the WavePropagation.visitParticlesWithId method to collect model
    particle data. Objects simply collect the visited particle
    data into a list.
    Used in conjunction with the
    esys.lsm.sim.WavePropagation.WavePropagation.visitParticlesWithId
    method.
    """
    def __init__(self):
        """
        Initialise with empty particle list.
        """
        self.particleList = []

    def __iter__(self):
        """
        Returns the iterator of the particle list.
        """
        return iter(self.particleList)

    def visitAParticle(self, particle):
        """
        Simply adds the specified particle data to the
        visited list.
        @type particle: L{esys.lsm.LsmPy.NRotSphere} or L{esys.lsm.LsmPy.RotSphere}
        @param particle: Particle data.
        """
        self.particleList.append(particle)

    def visitNRotSphere(self, p):
        self.visitAParticle(p)

    def visitRotSphere(self, p):
        self.visitAParticle(p)

class WavePropagation(WavePropagationPy.WavePropagation):
    """
    Wave propagation model, extends the
    esys.lsm.sim.WavePropagationPy.WavePropagation by providing
    convenient methods for saving data, creating partices/bonds and
    creating seismographs.
    """
    def __init__(
        self,
        domainBox,
        numWorkerProcesses = 2,
        mpiDimList = [0,0,0],
        do2d = False,
        timeStepSize = 0.05,
        dampingViscosity = 0.0
    ):
        """
        Initialise a wave propagation model.
        @type domainBox: L{BoundingBox}
        @param domainBox: a BoundingBox which specifies the model
                          particle domain.
        @type numWorkerProcesses: int
        @param numWorkerProcesses: the number of MPI worker processes used
                                   in model computations.
        @type mpiDimList: list
        @param mpiDimList: sequence of 3 int elements specifying the domain
                           decomposition used for the MPI worker processes.
        @type do2d: bool
        @param do2d: if True, forces the model to perform 2d calculations
        (particles move only in the x-y plane).
        @type timeStepSize: float
        @param timeStepSize: size of the time-step in the explicit integration
                             scheme.
        @type dampingViscosity: float
        @param dampingViscosity: if > 0 a damping viscosity is applied to
                                 particles.
        """
        WavePropagationPy.WavePropagation.__init__(
            self,
            numWorkerProcesses,
            mpiDimList
        )
        self.initVerletModel(
            "NRotSphere",
            2.5,
            0.20
        )
        self.setTimeStepSize(dt=timeStepSize)
        self.setSpatialDomain(domainBox)
        self.force2dComputations(do2d)
        if (dampingViscosity > 0.0):
            self.createDamping(
                DampingPrms(
                    type="Damping",
                    name="dampie-damp",
                    viscosity = dampingViscosity,
                    maxIterations=100
                )
            )
        self.sourceList      = []
        self.seismoGroups = SeismographGroupCollection(self)

    def createParticles(self, particles):
        """
        Creates model particles.
        @type particles: sequence
        @param particles: A sequence of L{esys.lsm.geometry.SimpleSphere}
                          objects.
        """
        WavePropagationPy.WavePropagation.createParticles(self, particles)

    def createBonds(self, connections, tagBondPrmDict):
        """
        Creates elastic bonds between particles.
        @type connections: sequence of
                           L{esys.lsm.geometry.GeometryPy.TaggedIdConnection}
        @param connections: create a bond for each connection which has
                            a tag key in the tagBondPrmDict dictionary.
        @type tagBondPrmDict: dict of tag:float
        @param tagBondPrmDict: a dictionary of (connectionTag:normalK)
                              pairs. Connections with tag connectionTag
                              will have a corresponding linear elastic
                              bond created with spring constant normalK.
        """
        self.createConnections(connections)
        for tag in list(tagBondPrmDict.keys()):
            if (isinstance(tagBondPrmDict[tag], float)):
                bondPrms = \
                    NRotBondPrms(
                        name = "BondageTag_" + str(tag),
                        normalK = tagBondPrmDict[tag],
                        breakDistance = 10.0,
                        tag = tag
                    )
            else:
                bondPrms = tagBondPrmDict[tag]
            self.createInteractionGroup(bondPrms)

    def addSource(self, source):
        """
        Add a source disturbance to the model.
        @type source: L{WaveSource}
        @param source: source disturbance object
        """
        self.sourceList.append(source)
        self.addPreTimeStepRunnable(source)

    def createSource(self, prm):
        """
        Creates a propagation source (a WaveSource object)
        within the model.
        @type prm: L{SourcePrms}
        @param prm: parameters for creating a WaveSource object.
        """
        self.addSource(WaveSource(prm, self))

    def getTime(self):
        """
        Returns simulation time value (simply number of time steps
        multiplied by the time step size).
        @rtype: float
        @return: Simulation time.
        """
        return self.getTimeStep()*self.getTimeStepSize()

    def createSources(self, sourcesPrms):
        """
        Creates multiple propagation sources.
        @type sourcesPrms: sequence of L{SourcePrms} objects
        @param sourcesPrms: A WaveSource is created for each object in
                            this sequence.
        """
        if (not hasattr(sourcesPrms, "__iter__")):
            sourcesPrms = [sourcesPrms]
        for prm in sourcesPrms:
            self.createSource(prm)

    def createSeismographGroup(
        self,
        posnIteratable,
        fileNamePrefix,
        sourcePosn
    ):
        """
        Creates a group of seismographs (a SeismographGroup obect)
        from a specified sequence of spatial locations.
        @type posnIteratable: sequence of L{Vec3}
        @param posnIteratable: a sequence of coordinates which indicate
                               the locations of seismographs. A Seismograph
                               object is created for each point in
                               posnIteratable.
        @type fileNamePrefix: string
        @param fileNamePrefix: Seismograph data is saved to files with this
                              file name prefix.
        @type sourcePosn: L{Vec3}
        @param sourcePosn: location of the point-source disturbance.
        """
        self.seismoGroups.createGroup(
            posnIteratable,
            fileNamePrefix,
            sourcePosn
        )

    def saveSeismoData(self):
        """
        Appends seismograph data to file for the current time step.
        """
        self.seismoGroups.saveSeismoData()

    def writeReorderedRecordSectionData(self):
        """
        Reads time-ordered record-section data from file
        and writes a distance-from-source ordered file.
        """
        self.seismoGroups.writeReorderedRecordSectionData()

    def getParticleDataFileName(self, idx, fileNamePrefix):
        """
        Returns the name of the file where particle data is saved
        for the specifed index.
        @type idx: int
        @param idx: The snap-shot/frame number
        @type fileNamePrefix: string
        @param fileNamePrefix: The returned file name is a concatenation
                              of this prefix with a suffix formed using
                              the index argument.
        """
        return fileNamePrefix + "{0:d}.txt".format(idx,)

    def writeParticleDataOrig(self, idList, index, fileNamePrefix="particle_"):
        """
        Pure python-implemented version of saving particle displacement
        data, the C++ implementation provided by
        esys.lsm.sim.WavePropagation.WavePropagation class is faster.
        Writes particle displacement and velocity data to file.
        Each line of the file is
        'px py pz dx dy dz vx vy vz' where p=(px,py,pz) is the particle
        position, d=(dx,dy,dz) is the current particle displacement
        (ie position relative to initial position) and (vx,vy,vz) is the
        particle velocity.
        @type idList: sequence of int
        @param idList: list of particle id's for which data will be saved.
        @type index: int
        @param index: integer used to generate file name.
        @type fileNamePrefix: string
        @param fileNamePrefix: prefix of file where displacement data is saved.
        """
        #
        # Collect data for all particles with id's specified
        # in the idList.
        #
        visitor = PVisitor()
        self.visitParticlesWithId(idList, visitor)
        f = file(self.getParticleDataFileName(index, fileNamePrefix), "w")
        for p in visitor:
            f.write(
                str(p.getPosn()) +\
                " " +\
                str(p.getPosn()-p.getInitialPosn()) +\
                " " +\
                str(p.getVelocity()) +\
                "\n"
            )
        f.close()

    def writeParticleData(self, index, fileNamePrefix="particle_"):
        """
        Writes particle displacement and velocity data to file.
        Each line of the file is
        'px py pz dx dy dz vx vy vz' where p=(px,py,pz) is the particle
        position, d=(dx,dy,dz) is the current particle displacement
        (ie position relative to initial position) and (vx,vy,vz) is the
        particle velocity. Particle data is saved for all particles with an
        id in the list returned by the self.getParticleDataIdList.
        @type index: int
        @param index: integer used to generate file name.
        @type fileNamePrefix: string
        @param fileNamePrefix: prefix of file where displacement data is saved.
        """
        #
        # Collect data for all particles with id's specified
        # in the idList.
        #
        self.writeParticleDataToFile(
            self.getParticleDataFileName(index, fileNamePrefix)
        )
