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
from esys.lsm.util     import Vec3, BoundingBox, seedDefaultRng
from esys.lsm.geometry import SimpleSphere,SimpleBlock, CubicBlock
from esys.lsm.geometry import HexagBlock, RandomBoxPacker, ParticleCollection
from esys.lsm.geometry import RandomSpherePacker

__initCfgTutSection = \
"""
Tutorial 2: Creating Initial Configurations of Particles
========================================================

Part of designating the initial conditions for a discrete element
simulation is to specify a starting configuration of particles.
The `esys.lsm.geometry` package provides classes for generating
basic initial *packings* of spherical particles. This tutorial
demonstrates the use of these classes for generating various
types of packings.

Regular Sphere Packings of Uniformly Sized Spheres
--------------------------------------------------

There are three classes for generating rectangular blocks of
regularly-packed particles:
  
  `esys.lsm.geometry.SimpleBlock`
    Spheres packed such that their centre points form
    a regular grid. Each sphere has a maximum of six
    contacting neighbours.
  `esys.lsm.geometry.CubicBlock`
    Spheres packed in a cubic close packing.
  `esys.lsm.geometry.HexagBlock`
    Spheres packed in a hexagonal close packing.

A good introduction to sphere packing, which defines
*cubic close* and *hexagonal close* packings, can be
found at http://mathworld.wolfram.com/SpherePacking.html.

The following Python code illustrates the generation of three
packing configurations::

  from esys.lsm.geometry import SimpleBlock,CubicBlock,HexagBlock

  simpleBlk = SimpleBlock(dimCount = [10,20,5], radius = 0.5)
  cubicBlk  = CubicBlock(dimCount = [20,20,20], radius = 0.25)
  hexagBlk  = HexagBlock(dimCount = [5,10,15], radius = 1.0)

The ``SimpleBlock``, ``CubicBlock`` and ``HexagBlock`` constructors
each take two arguments:
  
  ``dimCount``
    A three-element list of integers indicating the number of particles
    in each axis direction.
  ``radius``
    The radius of each particle in the packing.

Thus, the ``simpleBlk`` packing, consists of ``10*20*5=1000`` particles
of radius ``0.5`` arranged in an axis-aligned block which is ``10`` particles
wide in the *x* direction, ``20`` particles high in the *y* direction and
``5`` particles deep in the *z* direction. Similarly, the ``cubicBlk``
packing consists of ``20*20*20=8000`` particles of radius ``0.25`` and
the ``hexagBlk`` packing consists of ``5*10*15=750`` particles of radius
``1.0``. The ``cubicBlk`` and ``hexagBlk`` packings are produced by arranging
2D layers of spheres on top of one another. The 2D layers themselves are
produced by arranging linear rows of spheres in a regular triangular packing.
By default, each row is aligned with the *x* axis and the 2D layers are
generated parallel to the *x*-*z* plane. This default arrangement is referred
to as an ``XZY`` *orientation* (``XZ`` layers packed in the ``Y`` direction).
The ``HexagBlock`` and ``CubicBlock`` constructors accept an additional
``orientation`` keyword argument, for specifying alternative orientations.
For example, the above assignment statements involving the ``cubicBlk``
and ``hexagBlk`` variables are equivalent to::

  from esys.lsm.geometry import SimpleBlock,CubicBlock,HexagBlock,Orientation

  cubicBlk  = CubicBlock(dimCount = [20,20,20], radius = 0.25, orientation = Orientation.XZY)
  hexagBlk  = HexagBlock(dimCount = [5,10,15], radius = 1.0, orientation = Orientation.XZY)

Alternative ``orientation`` values produce different layer packings. For example::

  hexagBlk  = HexagBlock(dimCount = [5,10,15], radius = 1.0, orientation = Orientation.YZX)

specifies an ``YZX`` orientation. This produces a packing where ``10``
linear-rows (each row consisting of ``5`` particles) parallel with the
*y*-axis, are stacked into a layer which is parallel with the *y*-*z*
plane. ``15`` of these *y*-*z* layers are then packed in the *x* direction
to form the block.

Packings of Randomly Sized Spheres
----------------------------------

The following script generates a packing of randomly sized spheres
within an axis-aligned box::

  from esys.lsm.geometry import RandomBoxPacker
  from esys.lsm.util     import Vec3, BoundingBox, seedDefaultRng

  minR = 0.2
  maxR = 1.0
  seedDefaultRng(123454321)
  rndPkr = \\
    RandomBoxPacker(
        minRadius = minR,
        maxRadius = maxR,
        cubicPackRadius = maxR,
        maxInsertFails = 32000,
        bBox = BoundingBox(Vec3(0,0,0), Vec3(10,20,5)),
        circDimList = [False,False,False],
        tolerance = 0.01*minR
    )
  rndPkr.generate()
  rndBlk = rndPkr.getSimpleSphereCollection()
  print "Number of Spheres = ", len(rndBlk)

These statements generate a packing of spheres with random radius
values ranging between 0.2 and 1.0. After initialising the
minimum radius (``minR``) and maximum radius (``maxR``) variables,
a `RandomBoxPacker` object is created and assigned to the ``rndPkr``
variable. The ``RandomBoxPacker`` constructor accepts the following
arguments:
  
  ``minRadius``
    Minimum radius of generated spheres.
  ``maxRadius``
    Maximum radius of generated spheres.
  ``cubicPackRadius``
    Initial randomly sized seed particles
    are generated at positions corresponding to a regular
    cubic-packed grid of spheres with radius ``cubicPackRadius``.
  ``maxInsertFails``
    Stopping criterion for terminating the random
    insertion algorithm. Algorithm terminates after ``maxInsertFails``
    number of attempts at fitting a randomly generated particle into
    the box.
  ``bBox``
    A box specifying the region into which particles are packed.
    A 2D packing can be generated by specifying a zero sized *z* dimension,
    eg ``bBox=BoundingBox(Vec3(1,1,0),Vec3(21,21,0))`` will generate a 2D
    packing in the *x*-*y* plane.
  ``circDimList``
    A list of 3 boolean values indicating which (if any)
    of the box dimensions is circular (note, only a single dimension may
    be circular).
  ``tolerance``
    Generated particles may overlap by no more than this amount.

The construction of the ``RandomBoxPacker`` object does not produce
the packing, the packing is generated by the call to the
`RandomBoxPacker.generate` method (ie the ``rndPkr.generate()``
statement). The collection of random-sized particles
is assigned to the ``rndBlk`` variable. In this example, with
the above parameters, the ``rndBlk`` collection contains
``len(rndBlk) == 3501`` `SimpleSphere` objects.

The packing algorithm (executed in the ``RandomBoxPacker.generate`` method)
has two main parts: the generation of an initial group of *seed* spheres
and random insertion of *fill-in* spheres. The seed sphere radii are
selected from a uniformly random distribution, and the seed sphere centre
positions cooincide with the centre points of spheres of radius
``cubicPackRadius`` arranged in a regular cubic close packing.
Fill-in spheres are then generated by randomly choosing locations
within the specifed ``bBox`` box and attempting to find a sphere which
touches neighbouring spheres (or touches neighbouring spheres and a side
of the box) and which has radius in lying in the interval
``[minRadius,maxRadius]``. If no sphere can be found which fits with
neighbouring spheres and the boundary then this is recorded as a failed
insertion attempt and another random location is selected. When there
have been ``maxInsertFails`` *consecutive* failed insertion attempts,
the algorithm stops and the ``RandomBoxPacker.generate`` method returns.

Useful Methods for Manipulating `SimpleSphere` Collections
----------------------------------------------------------

The `SimpleBlock`, `CubicBlock` and `HexagBlock` are all subclasses
of the `ParticleCollection` class while the `RandomBoxPacker` and
`RandomSpherePacker` classes provide a ``getParticleCollection``
method for obtaining a ``ParticleCollection`` object.
An object of class ``ParticleCollection`` is simply a
list/group/collection of `SimpleSphere` objects. The ``ParticleCollection``
class provides methods which are useful for manipulating
the group of `SimpleSphere` objects. This section provides
some examples of using these methods.

Geometrical Transformations: Translation and Rotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once a basic arrangement of spheres has been generated it is
often useful to be able to translate, or rotate all of the
spheres in the collection. Objects of class `ParticleCollection`
provide the `ParticleCollection.translate` and
`ParticleCollection.rotate` methods which translate and rotate
all ``SimpleSphere`` objects contained in the collection::

  from __future__        import division
  from esys.lsm.geometry import SimpleBlock
  from esys.lsm.util     import Vec3
  from math              import pi
  
  smplBlk = SimpleBlock(dimCount=[10,3,1], radius=0.5)
  smplBlk.translate(translation = Vec3(-5, 0, 0))
  smplBlk.rotate(axis = Vec3(0, 0, pi/4.0), pt = Vec3(0, 0, 0))

Here the ``smplBlk`` variable is assigned to a collection of 30
``SimpleSphere`` objects arranged in a regular grid. The ``translate``
method translates all particles in the collection by ``5`` units in
the negative *x* direction. All the ``SimpleSphere`` objects are then
rotated by ``pi/4.0`` radians counterclockwise about the *z* axis.

Iterating
^^^^^^^^^

The ``ParticleCollection`` class supports the python *iterable*
concept, meaning that individual elements of the
collection can be enumerated in a loop as follows::
  
  from __future__        import division
  from esys.lsm.geometry import SimpleBlock
  from esys.lsm.util     import Vec3
  from math              import pi
  
  smplBlk = SimpleBlock(dimCount=[10,3,1], radius=0.5)
  for ss in smplBlk:
      ss.translate(Vec3(-5,0,0))
      ss.rotate(axis = Vec3(0, 0, pi/4.0), pt = Vec3(0, 0, 0))

During each iteration in the ``for`` loop, the *loop variable*
``ss`` is assigned to a ``SimpleSphere`` object in the ``smplBlk``
collection. The above statements result in a final state which is
equivalent to the transformations performed by directly calling
methods of the ``smplBlk`` object. That is, the above ``for``
loop is equivalent to::
  
  smplBlk = SimpleBlock(dimCount=[10,3,1], radius=0.5)
  smplBlk.translate(translation = Vec3(-5, 0, 0))
  smplBlk.rotate(axis = Vec3(0, 0, pi/4.0), pt = Vec3(0, 0, 0))

Saving to File
^^^^^^^^^^^^^^

Often it is the case that it is desired to run multiple simulations
with the same initial configuration of spheres. Rather than generating
the packing for each simulation, it is convenient to save the initial
configuration of spheres to file and read in this configuration
before executing the simulation. Saving the configuration to file is
also convenient when it is time consuming to generate the initial
configuration, that is, often it is much faster to read a packing
of spheres from file rather than compute the positions of particles
in the packing each time the data is needed.

There are two options for *serializing* the data: custom file formats
and Python's built-in serialization modules
pickle_ and cPickle_.

.. _pickle: http://www.python.org/doc/2.4.2/lib/module-pickle.html
.. _cPickle: http://www.python.org/doc/2.4.2/lib/module-cPickle.html

Custom File Formats
~~~~~~~~~~~~~~~~~~~

Python's I/O interfaces and simple conversion from basic types
to strings (and vice versa) make it relatively simple to save
any kind of data to file::

  from esys.lsm.geometry import HexagBlock
  
  hexagBlk = HexagBlock(dimCount=[10,10,10], radius=0.25)
  f = open("HexagBlk.txt", "w")
  for ss in hexagBlk:
      f.write(
          "{0!s} {1!s} {2!s} {3!s}\\n"\\
          .format(ss.getId(), ss.getRadius(), ss.getCentre(), ss.getMass())
      )
  f.close()

After constructing the hexagonal close packing of 1000 spheres
(and assigning it to the ``hexagBlk`` variable),
the file_ ``"HexagBlk.txt"`` is opened in *write* mode.
The ``for`` loop then iterates over each ``SimpleSphere`` object
in the ``hexagBlk`` collection and writes a line to the
``"HexagBlk.txt"`` file. The ``"{0!s} {1!s} {2!s} {3!s}\\n".format(...)``
argument to the ``f.write`` method call is a `string formatting
operation`.

__ http://www.python.org/doc/2.4.2/lib/typesseq-strings.html
.. _file: http://www.python.org/doc/2.4.2/lib/built-in-funcs.html#l2h-25

The data can be read from file to reconstruct a ``HexagBlock``
as follows::

  from esys.lsm.geometry import HexagBlock, SimpleSphere
  from esys.lsm.util     import Vec3
  
  hexagBlk = HexagBlock()
  for line in file("HexagBlk.txt", "r"):
      elemList = line.split(" ")
      hexagBlk.createParticle(
        SimpleSphere(
            id     = int(elemList[0]),
            radius = float(elemList[1]),
            centre = Vec3(map(float,elemList[2:5])),
            mass   = float(elemList[5])
        )
      )

First the ``hexagBlk`` variable is assigned to an empty
``ParticleCollection`` object. The ``for`` loop iterates over
each line in the ``"HexagBlk,txt"`` file. The ``line`` string
is split into a list of string elements, which are delimited
by the ``" "`` space character, and this list is assigned to
the ``elemList`` variable. In the next statement, a ``SimpleSphere``
object is created within the ``hexagBlk`` collection by calling
the ``hexagBlk.createParticle`` method. The arguments to the
``SimpleSphere`` constructor are converted from strings to either
an ``int`` value or ``float`` values. The call to
``map(float, elemList[2:5])`` converts the ``elemList[2:5]``
list-slice of three string elements to a list of three ``float``
elements.

Saving data in a custom file format offers the advantage of
easily being able to observe the data in a text file. Explicity
saving data in a specific format may also make it possible to use
other software to visualise/analyse the data.

Pickling
~~~~~~~~

Python's pickle_ and cPickle_ modules provide a convenient mechanism
for serializing Python objects. In Python-speak, serialization
is called *pickling* and de-serialization is called *unpickling*.
The following Python code pickles the same data as that of the
`Custom File Formats`_ section::

  from esys.lsm.geometry import HexagBlock
  import pickle
  
  hexagBlk = HexagBlock(dimCount=[10,10,10], radius=0.25)
  f = file("HexagBlk.pkl", "w")
  pickle.dump(obj = hexagBlk, file = f, protocol = pickle.HIGHEST_PROTOCOL)
  f.close()

As usual, the ``hexagBlk`` is assigned to a hexagonal close-packed
collection of ``SimpleSphere`` objects. The file ``HexagBlk.pkl``
is opened in *write* mode, and the call to the ``pickle.dump`` function
takes care of pickling the ``hexagBlk`` object to file.
The data is unpickled by executing::
  
  import pickle
  
  f = file("HexagBlk.pkl", "r")
  hexagBlk = pickle.load(file = f)
  f.close()

The ``HexagBlk.pkl`` file is opened in *read* mode and
the call to ``pickle.load`` unpickles the data from
the file into the original pickled object and this
is assigned to the ``hexagBlk`` variable.

The advantage of the pickling technique is that it requires far
fewer lines of code and the user doesn't have to know the internal
structure of the data being saved. The main disadvantage is that
the saved data-file is not human-readable and is only loadable
by Python.

"""


__doc__ = \
    esys.lsm.doc.Util.setSectionDoc("InitialConfigSection",__initCfgTutSection) \
    + "\n:summary: Creating packings of spheres and saving/loading packing-data from file.\n"
