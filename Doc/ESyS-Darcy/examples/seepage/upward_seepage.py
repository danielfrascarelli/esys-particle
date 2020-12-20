#upward_seepage.py: A upward seepage simulation through particle assembly using ESyS-Darcy (a fluid-coupled version of ESyS-Particle)
#       Author: Q.Shao
#       Date:  21 August 2017
#       Organisation: University of Queensland
#       (C) All rights reserved, 2017.
#
#
				
#import the appropriate ESyS-Particle modules:
from esys.lsm import *
from esys.lsm.util import *
from esys.lsm.geometry import *

#instantiate a simulation object
#and initialise the neighbour search algorithm:
sim = LsmMpi(numWorkerProcesses=1, mpiDimList=[1,1,1])

sim.initNeighbourSearch(
    particleType="RotSphere",
    gridSpacing=0.008,
    verletDist=0.0003
)

#specify the number of timesteps and timestep increment:
sim.setNumTimeSteps(500000)
sim.setTimeStepSize(1.0e-6)

#specify the spatial domain for the simulation:
domain = BoundingBox(Vec3(0.,0.,-0.01), Vec3(0.04,0.02,0.07))
sim.setSpatialDomain(domain)

#add fluid to the domain:
sim.createFluidInteraction(
  cellside=0.01,         #side length of fluid cell
  Bw=1.0e6,              #bulk modulus of fluid
  Bp=1.0e6,              #bulk modulus of particle
  Mu=1.0e-3,             #fluid viscosity
  alpha=0.5,             #adjusting factor between two time steps
  flowrate=0.1,          #rate of inflow
  pressure=0.,           #hydraulic pressure
  inflow=Vec3(0,0,1),    #directions of inflows
  outflow=Vec3(0,0,1),   #directions of outflows 
  fluid_timestep=1.0e-4  #time step size for updating fluid phase
)

#read in geometry input file
sim.readGeometry("input.geo")

#set partilce density
sim.setParticleDensity (
  tag = 0,
  mask = -1,
  Density = 2600
)

#add a left side wall to the model:
sim.createWall (
  name = "left_wall",
  posn = Vec3(0.0000, 0.0000, 0.0000),
  normal = Vec3(1.0000, 0.0000, 0.0000)
)

#add a right side wall to the model:
sim.createWall (
  name = "right_wall",
  posn = Vec3(0.04, 0.0000, 0.0000),
  normal = Vec3(-1.0000, 0.0000, 0.0000)
)

#add a back side wall to the model:
sim.createWall (
  name = "back_wall",
  posn = Vec3(0.0000, 0.0000, 0.0000),
  normal = Vec3(0.0000, 1.0000, 0.0000)
)

#add a front side wall to the model:
sim.createWall (
  name = "front_wall",
  posn = Vec3(0.0000, 0.02, 0.0000),
  normal = Vec3(0.0000, -1.0000, 0.0000)
)

#add a bottom side wall to the model:
sim.createWall (
  name = "bottom_wall",
  posn = Vec3(0.0000, 0.0000, 0.0000),
  normal = Vec3(0.0000, 0.0000, 1.0000)
)

#specify that particles undergo frictional interactions:
sim.createInteractionGroup (
  RotFrictionPrms (
    name = "friction",
    normalK = 1.0e6,
    dynamicMu = 0.4,
    staticMu = 0.6,
    shearK = 1.0e6,
    scaling = True,
    rigid = True,
    meanR_scaling = True
  )
)

#specify that particles undergo elastic repulsion
#from the left side wall:
sim.createInteractionGroup (
  NRotElasticWallPrms (
    name = "lw_repel",
    wallName = "left_wall",
    normalK = 1.0e6
  )
)

#specify that particles undergo elastic repulsion
#from the right side wall:
sim.createInteractionGroup (
  NRotElasticWallPrms (
    name = "rw_repel",
    wallName = "right_wall",
    normalK = 1.0e6
  )
)

#specify that particles undergo elastic repulsion
#from the back side wall:
sim.createInteractionGroup (
  NRotElasticWallPrms (
    name = "bkw_repel",
    wallName = "back_wall",
    normalK = 1.0e6
  )
)


#specify that particles undergo elastic repulsion
#from the front side wall:
sim.createInteractionGroup (
  NRotElasticWallPrms (
    name = "fw_repel",
    wallName = "front_wall",
    normalK = 1.0e6
  )
)


#specify that particles undergo elastic repulsion
#from the bottom side wall:
sim.createInteractionGroup (
  NRotElasticWallPrms (
    name = "btw_repel",
    wallName = "bottom_wall",
    normalK = 1.0e6
  )
)

#add a CheckPointer to store simulation data
sim.createCheckPointer (
  CheckPointPrms (
  fileNamePrefix = "results/snap",
    beginTimeStep=0,
    endTimeStep=500000,
    timeStepIncr=100000
  )
)

#create a FieldSaver to store distributed pore pressures
sim.createFieldSaver (
  FluidScalarFieldSaverPrms(
    fieldName="disP",
    fileName="results/Pressure",
    fileFormat="VTI",
    beginTimeStep=0,
    endTimeStep=500000,
    timeStepIncr=100000
  )
)

#create a FieldSaver to store porosities
sim.createFieldSaver (
  FluidScalarFieldSaverPrms(
    fieldName="Phi",
    fileName="results/Porosity",
    fileFormat="VTI",
    beginTimeStep=0,
    endTimeStep=500000,
    timeStepIncr=100000
  )
)

#create a FieldSaver to store permeabilities
sim.createFieldSaver (
  FluidScalarFieldSaverPrms(
    fieldName="K",
    fileName="results/Permeability",
    fileFormat="VTI",
    beginTimeStep=0,
    endTimeStep=500000,
    timeStepIncr=100000
  )
)

#Execute the simulation:
sim.run()
