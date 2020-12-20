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
__docformat__ = "restructuredtext en"

import esys.lsm.doc.Util
from esys.lsm.util import Vec3,BoundingBox
from esys.lsm      import LsmMpi, NRotSphere, GravityPrms

__gravTutSection = \
"""
Tutorial 1: A Particle Falling Due To Gravity
==============================================

This tutorial analyses a simple model where a particle falls
due to the influence of gravity. There are no fixed objects,
and no interactions between particles, just a spherical
particle falling in space. Although the example is very simple,
it does illustrate the important basic features of running a
simulation.

The Complete Script
-------------------

The complete script for the simulation is as follows::

  from esys.lsm      import *
  from esys.lsm.util import Vec3, BoundingBox
  
  sim = LsmMpi(numWorkerProcesses=4, mpiDimList=[2,2,1])
  sim.initVerletModel(
      particleType = "NRotSphere",
      gridSpacing  = 2.5,
      verletDist   = 0.5
  )
  
  domain = BoundingBox(Vec3(-20,-20,-20), Vec3(20, 20, 20))
  sim.setSpatialDomain(domain)
  
  particle = NRotSphere(id=0, posn=Vec3(0,0,0), radius=1.0, mass=1.0)
  sim.createParticle(particle)
  
  sim.createInteractionGroup(
      GravityPrms(name="earth-gravity", acceleration=Vec3(0,-9.81,0))
  )
  
  sim.setTimeStepSize(0.001)
  sim.setNumTimeSteps(10000)
  sim.run()

Before explaining the above statements, the command line for running
this script is demonstrated.

Executing an LSM Simulation Python Script
-----------------------------------------

If the above script is saved to the file ``GravitySim.py``, the linux/unix
command line for running the script under LAM MPI is::

  $ mpiexec -machinefile hosts.txt -np 1 /absolute/path/to/lsm/mpipython GravitySim.py

where the file ``hosts.txt`` contains information about the various MPI
hosts. In its simplest form, the ``hosts.txt`` file can contain just a
single line of text::

  localhost

which causes all MPI processes to be invoked on the local machine.

To execute the script, e.g. on an SGI Altix, using SGI's MPI implementation
(MPT) the command line is::

  $ mpirun -up 5 -np 1 /absolute/path/to/lsm/mpipython GravitySim.py

In both these cases the initial number of MPI processes is specified
by the ``-np`` option and should always have value ``1``. The python
script is used to specifiy a number of dynamically created MPI processes
and the simulation computations are distributed across these dynamically
created processes. The MPT ``mpirun`` command requires the ``-up`` option
to specify the *MPI universe size*, which is just the number of all
MPI processes.

Analysing the Script
--------------------

Importing Packages and Modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to use the LSM python library to run simulations, it must be
*imported* into the script. The first two lines::
  
  from esys.lsm      import *
  from esys.lsm.util import Vec3, BoundingBox

achieve the import. The first statement imports a multitude
of classes and functions into the current python *namespace*
from the `esys.lsm` *package*.
The second ``import`` statement imports the two specified classes
into the the namespace: `Vec3` and `BoundingBox`. Objects of class
``Vec3`` are 3-element vectors of floating point values (representing
3D points, for example) and objects of class ``BoundingBox``
represent an axis-aligned rectangular-prism/box.

There are a number of subpackages defined within the `esys.lsm`
top-level package. Each subpackage defines classes and functions
which may be imported to provide various functionality:

  `esys.lsm.benchmarks`
    Defines modules and packages which in turn define simlation scripts
    which are used for performance analysis.
  
  `esys.lsm.doc`
    Defines documentation-only python modules with scripting tutorials.
  
  `esys.lsm.examples`
    Defines sub-packages and modules of example simulation scripts.
  
  `esys.lsm.geometry`
    Defines modules for generating initial packings/configurations of
    LSM particles and bonds.

  `esys.lsm.sim`
    Base classes for different types of simulations and helper classes
    for running *suites* of simulations.

  `esys.lsm.util`
    Utility package providing basic common functionality.

  `esys.lsm.vis`
    Visualisation packages and modules.

Creating and Initialising the Simulation Object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the LSM modules have been imported, the simulation object
is created and initialised::

  sim = LsmMpi(numWorkerProcesses=4, mpiDimList=[2,2,1])
  sim.initVerletModel(
      particleType = "NRotSphere",
      gridSpacing  = 2.5,
      verletDist   = 0.5
  )

The first statement creates an `LsmMpi` object. ``LsmMpi`` objects
provide the means to define and run LSM simulations. Here, two
arguments are provided to the ``LsmMpi`` *constructor*:

  ``numWorkerProcesses``
    The number of MPI *worker* processes which are
    dynamically created.

  ``mpiDimList``
    A 3 element list of integers which specify a regular
    axis-aligned grid decomposition of the rectangular domain. Each
    domain-cell of the grid is assigned to an MPI process. The individual
    MPI processes are then responsible for performing calculations on
    particles which lie in their assigned cell.

The second statement, is a call to the `LsmMpi.initVerletModel` *method*
of the ``sim`` object. This method defines the type of discrete-element
particles which are to be used in the simulation and also initialises
data-structures used in the contact detection algorithm. The three arguments
are defined as follows:

  ``particleType``
    A string defining the type of discrete-element particle. In this example,
    particles are non-rotational spheres ``"NRotSphere"``.
  
  ``gridSpacing``
    A regular grid of cubic cells is used to efficiently determine the
    neighbouring particles of each particle. This argument defines the
    length of the cubic cell sides and needs to be greater than the
    maximum particle radius.
  
  ``verletDist``
    This distance determines the frequency with which neighbour-lists/contacts
    are updated. If any particle moves further than this distance during the
    course of a simulation, then the neighbour lists (and interactions/contacts)
    are updated.

Optimal values for ``gridSpacing`` and ``verletDist`` are simulation
dependent and the parameters must obey the constraint
``gridSpacing > maxRadius + verletDist``, where ``maxRadius`` is
the maximum radius of all particles. Smaller ``verletDist`` values
will lead to smaller neighbour/contact lists and faster force calculations.
However, a small ``verletDist`` value also means more frequent
recalculation of the neighbour lists.


Setting the Spatial Domain
^^^^^^^^^^^^^^^^^^^^^^^^^^

The next step in parameterising the simulation is to specify the spatial
domain over which particles are *tracked*. The LSM allows the specification
of a rectangular box as the bounding domain::
  
  domain = BoundingBox(Vec3(-20,-20,-20), Vec3(20, 20, 20))
  sim.setSpatialDomain(domain)

Here, the ``domain`` variable has been set as a `BoundingBox`
object. The ``BoundingBox`` constructor takes two `Vec3` arguments,
a minimum point (lower left back corner) ``Vec3(-20,-20,-20)``
and a maximum point (upper right front corner) ``Vec3(20, 20, 20)``.
This domain specifies a cube with side length ``40`` and centred at
the origin.

The `LsmMpi.setSpatialDomain` method is called to specify the
simulation domain. Equivently, the domain could have been set
in a single statement as follows::
  
  sim.setSpatialDomain(BoundingBox(Vec3(-20,-20,-20), Vec3(20, 20, 20)))
    
Creating a Spherical Particle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now that the internal contact detection object and spatial domain have
been initialised, objects relating to the actual discrete element
physics can be created. First a particle is constructed within the model::

  particle = NRotSphere(id=0, posn=Vec3(0,0,0), radius=1.0, mass=1.0)
  sim.createParticle(particle)

The ``particle`` variable is assigned to a `NRotSphere` object and
this object is then provided as an argument to the `LsmMpi.createParticle`
method. Inside the ``sim.createParticle`` method, an identical **copy**
of the ``particle`` argument is constructed and added internally to the
``LsmMpi`` simulation ``sim`` object.

Creating the Gravity Interaction Group
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to create simulations of interest, not only must the
initial particle conditions of a simulation be specified but also,
forces imparted on particles also should to be introduced.
Within the LSM scripting environment, these different types of forces
are termed *interaction groups*. Creating an interaction group within
the model specifies how discrete-element particles interact with each
other and also how they interact with other entities within the model.
Gravity is introduced into the model by creating an interaction group
as follows::
  
  sim.createInteractionGroup(
      GravityPrms(name="earth-gravity", acceleration=Vec3(0,-9.81,0))
  )

The `LsmMpi.createInteractionGroup` method, of the ``sim`` object, is called
with a `GravityPrms` argument. The ``GravityPrms`` argument specifies
that a gravity interaction group is created in the model, and this gravity
group causes all particles in the model to be subjected to a *gravitational*
force with acceleration ``acceleration=Vec3(0,-9.81,0)``. All interaction
groups have an associated `InteractionPrms` subclass, and the constructors of
these ``InteractionPrms`` subclasses always have a ``name`` argument.
This ``name`` argument is used to uniquely identify an interaction group
when manipulating the group using the ``LsmMpi`` interface.

Running the Simulation
^^^^^^^^^^^^^^^^^^^^^^

To execute the time-stepping scheme for integrating the equations
of motion, the time-step size and number of time-steps are
specified before calling the `LsmMpi.run` method of the ``sim``
object::

  sim.setTimeStepSize(0.000025)
  sim.setNumTimeSteps(600000)
  sim.run()

Instead of executing all time-steps via the ``run`` method,
individual steps may be specified in a loop as follows::

  sim.setTimeStepSize(0.000025)
  maxTime = 600000*sim.getTimeStepSize()
  t = 0.0
  while (t < maxTime):
      sim.runTimeStep()
      t += sim.getTimeStepSize()

Here, the `LsmMpi.runTimeStep` method is called within a
``while`` loop to execute individual time-steps. This form
of performing the integration allows the execution of arbitrary
additional code before and after each time-step.

Summary
-------

This simple gravity simulation example illustrates the basic steps
involved in all LSM simulations scripts:

  1. Construct and initialise an LSM simulation object (`LsmMpi`).
  2. Specify the spatial domain.
  3. Create an initial configuration of particles within the model.
  4. Specify inter-particle interactions and/or interactions between
     particles and the environment.
  5. Executing the time integration.

Of course, the purpose of running simulations is to produce data.
In this simple gravity example, neither saving data or visualisation
of particles has been covered. These topics are examined in subsequent
tutorials.

"""

__doc__ = \
    esys.lsm.doc.Util.setSectionDoc("GravityTutSection",__gravTutSection) \
    + "\n:summary: A simple gravity simulation.\n"
