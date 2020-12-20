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


from __future__ import division
import esys.lsm.doc.Util
from esys.lsm.vis.core import Scene, Camera, Sphere, GlyphData
from esys.lsm.geometry import SimpleSphere, SimpleBlock, CubicBlock
from esys.lsm.geometry import HexagBlock, RandomBoxPacker, ParticleCollection
from esys.lsm.geometry import RandomSpherePacker

__visualisationTutSection = \
"""
Tutorial 3: Visualisation
=========================

The `esys.lsm.vis` package provides sub-packages and modules for
rendering data to screen and/or image files. This tutorial covers
the fundamentals of the API.

POV-Ray and VTK Renderers
-------------------------

The `esys.lsm.vis` package provides a Python API which
delegates calls to 3rd-party *render* software.
Currently, there are two *renderer* modules, ``esys.lsm.vis.povray``,
which delegates calls to POV-ray_ ray-tracing software, and
``esys.lsm.vis.vtk`` which delegates calls to Kitwares VTK_ scientific
visualisation library.

.. _POV-Ray: http://www.povray.org
.. _VTK: http://www.vtk.org

The main reason for the existence of the ``esys.lsm.vis`` package
is to allow the same Python code to be executed to produce images
using either VTK or POV-Ray as the *render engine*.
The ``esys.lsm.vis`` API also hopefully simplifies visualisation
tasks with a minimal loss in flexibility.

Scenes, Cameras and Objects
---------------------------

There are four basic steps when creating a visualisation of a data-set:
  
  (1) Creating a `Scene` object
  (2) Adding objects to the Scene
  (3) Positioning the `Camera` to view the objects
  (4) Rendering the ``Scene`` to produce an image

Here is an example script which renders two spheres::

  import esys.lsm.vis.povray
  import esys.lsm.vis.vtk

  renpkg = esys.lsm.vis.povray

  scene = renpkg.Scene()
  scene.add(renpkg.Sphere(center=[0.0, 0.0, 0.0], radius=1.0))
  scene.add(renpkg.Sphere(center=[2.0, 0.0, 0.0], radius=0.75))
  scene.getCamera().setPosn([0.0, 0.0, -4.0])
  scene.getCamera().setLookAt([1.0, 0.0, 0.0])
  scene.render(size=[512,384], offScreen=False, interactive=True)
  
As usual, the script begins with ``import`` statements to load
packages and define classes/functions. Here, the
``esys.lsm.vis.povray`` and ``esys.lsm.vis.vtk``
modules are both imported into the script. Next, the
``esys.lsm.vis.povray`` package is assigned to the ``renpkg``
variable and the classes from this package will be
used to construct the scene and generate an image.
An object of class ``renpkg.Scene`` is created and
assigned to the ``scene`` variable. ``Scene`` objects
are containers for objects appearing in a visualisation.
Two calls to the `Scene.add` method are then used to add two
`Sphere` objects to the scene. The camera is then
positioned at coordinate ``(0,0,-4)`` with a call
to the `Camera.getPosn` method and the camera is *pointed at*
the coordinate ``(1,0,0)`` which is situated between the
two spheres. Lastly, the scene is rendered as an image
with the call to the `Scene.render` method. The call to
the ``render`` method is passed three arguments:
  
  ``size``
    A pair of ``int`` values specifying the size of the image.
    In the above example, the image size is 512 pixels wide and
    384 pixels high.
  ``offScreen``
    A ``bool`` indicating whether a window is displayed to screen.
    In the above example, this argument is set to ``False``
    so that a window is displayed on-screen.
  ``interactive``
    A ``bool`` indicating whether the displayed window allows
    interaction. Typically, this means that the ``render`` method
    will *pause* until the user exits the window (by typing the
    letter ``q``). The ``esys.lsm.vis.vtk.Scene`` renderer window
    allows the user to manipulate the camera (zoom, translate and rotate)
    with the mouse. The ``interactive`` argument only has effect if
    ``offScreen`` is set to ``True``.

The ``Scene.render`` method also accepts ``fileName`` and ``imageType``
arguments for saving an image to file, see `Scene.render` documentation
for details.

Simple Animations
-----------------

Here, the camera is rotated about the look-at point
to generate a simple animation::
  
  import esys.lsm.vis.povray
  import esys.lsm.vis.vtk
  import math
  
  renpkg = esys.lsm.vis.vtk
  
  scene = renpkg.Scene()
  scene.add(renpkg.Sphere(center=[-3.5, 0.0, 0.0], radius=1.75))
  scene.add(renpkg.Cube(minPt=[-1.0, 0.0, 0.0], sideLength=2.75))
  scene.getCamera().setPosn([2.0, 6.0, -6.0])
  scene.getCamera().setLookAt([0.0, 0.0, 0.0])
  numFrames = 300
  deltaRadians = 2.0*math.pi/numFrames
  for i in xrange(0,numFrames):
      scene.getCamera().rotatePosn(
        axis = [0, deltaRadians, 0],
        axisPt = scene.getCamera().getLookAt()
      )
      scene.render(size=[512,384], offScreen=False, interactive=False)
  
  scene.render(size=[512,384], offScreen=False, interactive=True)

After the imports, VTK is chosen as the renderer package
(that is, ``esys.lsm.vis.vtk`` assigned to ``renpkg``),
and two objects are added to the scene: a sphere and a cube.
The initial position and look-at position of the camera is set
before entering the ``for`` loop to generate the frames. The
`Camera.rotatePosn` method is used to incrementally rotate the
camera position about a *y*-vector which passes through the camera
look-at coordinate. The scene is rendered on-screen in each step of
the ``for`` loop with the ``interactive`` parameter set to ``False``.
At completion of the ``for`` loop, the final frame is rendered once
again, but this time with ``interactive=True`` so that the script
pauses until the render window is exited.

Visualising a Data-Set Using GlyphData
--------------------------------------

An alternative to creating individual objects within a scene
as was done with the ``scene.add(Sphere(centre=(0,0,0), radius=1.0))``
statements is to use a `GlyphData` object. ``GlyphData`` objects
are used to generate multiple *glyphs*  (spheres, arrows, cylinders,
disks,...) of the same type for a given set of data. The following
example script uses a ``GlyphData`` object to display disks with
orientation information from a table of values::
  
  import esys.lsm.vis.povray
  import esys.lsm.vis.vtk

  #
  # Create table of data
  #
  #
  table = [
      # x, y, z,  radius, degrees
      [ 0, 0, 0,    0.20,    10.0]
      [ 1, 0, 0,    0.30,    40.0]
      [ 1, 1, 0,    0.40,    60.0]
      [ 0, 1, 0,    0.50,    80.0]
      [ 2, 0, 0,    0.40,   100.0]
      [ 2, 1, 0,    0.30,   120.0]
      [ 2, 2, 0,    0.30,   140.0]
  ]
  numRows = len(table)

  renpkg = esys.lsm.vis.povray
  diskExtractor = \\
      renpkg.DiskExtractor(
        radiusMap   = lambda idx: table[idx][3],
        centerMap   = lambda idx: table[idx][0:3],
        radiusScale = 1.0
      )
  disks = renpkg.GlyphData(xrange(0, numRows), diskExtractor)

  arrowExtractor = \\
      renpkg.ArrowExtractor(
          tailPtMap = lambda idx: table[idx][0:3]
          vecMap    = lambda idx: getVector(table[idx][3], table[idx][0:3], table[idx][4])
      )
  arrows = renpkg.GlyphData(xrange(0, numRows), arrowExtractor)

  scene = renpgk.Scene()
  scene.add(disks)
  scene.add(arrows)
  scene.getCamera().setPosn((1,1,-5))
  scene.getCamera().setLookAt((1,1,0))
  scene.render(offScreen=False, interactive=True, size=[512,512])

"""

__doc__ = \
    esys.lsm.doc.Util.setSectionDoc("VisualisationSection",__visualisationTutSection) \
    + "\n:summary: Visualising data tutorial.\n"
