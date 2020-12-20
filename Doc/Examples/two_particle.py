# import the esysparticle modules
from esys.lsm import *
from esys.lsm.util import *
from math import *

Rmin = 0.00442     
Rmax = 0.05       

#Create a simulation container:
sim = LsmMpi (
   numWorkerProcesses = 1,
   mpiDimList = [1,1,1]
)

#Initialise the neighbour search algorithm:
sim.initNeighbourSearch (
   particleType = "NRotSphere",
   gridSpacing = 2.5*Rmax,
   verletDist = 0.1*Rmin
)

#Set the number of simulation timesteps and the timestep size:
sim.setNumTimeSteps (1000)
sim.setTimeStepSize (0.0001)

#Define the spatial domain boundaries (min/max points):
domain = BoundingBox (
   minPt = Vec3 (-2,-2,-2),
   maxPt = Vec3 (2,2,2)  
)
sim.setSpatialDomain(domain)

density = 3000.
sep = Rmax

# create the first particle:
r1 = Rmax
m1 = 4.*3.141592654*density*r1**3./3.
particle1 = NRotSphere (
   id = 0,
   posn = Vec3 (-1.0*r1-0.5*sep, 0, 0),
   radius = r1,
   mass = m1
)
#Set the initial velocity of particle one:
particle1.setLinearVelocity ( Vec3 (0.5,0,0) )
particle1.setTag(1234)

#Insert the particle into the simulation container:
sim.createParticle (particle1)

#Create and insert a second particle:
r2 = Rmin
m2 = 4.*3.141592654*density*r2**3./3.
particle2 = NRotSphere (id=1, posn=Vec3(r2+0.5*sep,0,0), radius=r2, mass=m2)
particle2.setLinearVelocity ( Vec3(-0.5,0,0) )
particle2.setTag(1)
sim.createParticle(particle2)

#Specify frictional linear elastic repulsion between touching particles:
sim.createInteractionGroup (
   NRotFrictionPrms (
      name = "elastic_repulsion",
      normalK = 100000.0,
      dynamicMu = 0.6,
      shearK = 100000.0,
      scaling = True
   )
)

#add dashpot damping between touching particles:
sim.createInteractionGroup (
   LinearDashpotPrms (
      name = "damping",
      damp = 9.0,
      cutoff = 1.0
   )
)

#Create a CheckPointer to store snapshots of simulation data for visualisation:
# (e.g. store snapshots every 100 timesteps between timestep numbers 0 and 3000)
sim.createCheckPointer (
   CheckPointPrms (
      fileNamePrefix = "snapshot",
      beginTimeStep = 0,
      endTimeStep = 3000,
      timeStepIncr = 100
   )
)

#Execute the simulation:
outfile = open("data.csv","w")
N_max = sim.getNumTimeSteps()
for n in range (N_max):
   sim.runTimeStep()

   particles = sim.getParticleList()
   outfile.write("%d %f %f %f %f %f %f\n" % (n,particles[0].getPosn()[0],particles[0].getPosn()[1],particles[0].getPosn()[2],particles[1].getPosn()[0],particles[1].getPosn()[1],particles[1].getPosn()[2]))

outfile.close()
sim.exit()
