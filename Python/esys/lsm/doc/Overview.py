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
import esys.lsm.util.InstallInfo

__esysParticleOverviewSection = \
"""
{pkgName:s} Overview
==================================

The Lattice Solid Model (LSM) [MoraPAGEOPH1994]_, [PlaceJCompPhys1999]_ is a particle-based
model similar to the Distinct Element Model [CundallGeotech1979]_.
The model consists of particles which are characterized by
their shape, mass, position, orientation and velocity. The particles
interact with their nearest neighbours by imparting contact forces.
Typically, Discrete Element Model (DEM) particles are spherical and the
contact forces consist of a linear elastic normal component and linear
elastic tangential component. {pkgName:s} is a parallel implementation
of the LSM with a Python_ scripting interface.

.. _Python: http://www.python.org

Particle Types
--------------

Currently, there exist three types of {pkgName:s} spherical particles:
*non-rotational*, *rotational*, *thermal-rotational*:

  Non-rotational Spheres
    Non-rotational spherical particles possess no rotational
    degrees of freedom. Objects of class `esys.lsm.NRotSphere`
    represent non-rotational spherical particles.

  Rotational Spheres
    Rotational spherical particles possess orientation information.
    Particles of this type change orientation according to the
    applied moments. Objects of class `esys.lsm.RotSphere`
    represent rotational spherical particles.

  Thermal Rotational Spheres
    Thermal rotational spherical particles are the same as "Rotational Spheres"
    with the addition of thermal properties (temperature and thermal expansion).
    Objects of class `esys.lsm.RotThermalSphere`
    represent thermal rotational spherical particles.

Inter-particle Interactions
---------------------------

Interactions between model particles are also classified as
*non-rotational* and *rotational*. Two spherical particles
involved in a non-rotational interaction have
all forces applied at the centre of mass.
Two spherical particles involved in a
rotational interaction experience moments
due to forces which are, typically, applied
at a *contact point*. The inter-particle interactions
include:

  Non-rotational Elastic
    Purely linear elastic repulsion of particles when in contact.

  Non-rotational Bonded
    Linear elastic attraction and repulsion while bond remains intact.
    Bond *breaks* when a threshold separation distance is reached.

  Non-rotational Friction
    Linear elastic repulsion, linear elastic shear rigidity and Coulomb
    dynamic friction law.

  Rotational Elastic
    Linear elastic repulsion as well as linear elastic shear rigidity.

  Rotational Bonded
    Linear elastic tension, compression, shear, torsion and bending
    forces while bond remains intact. Bond *breaks* if a threshold
    force limit is reached.

  Rotational Friction
    Linear elastic repulsion, linear elastic shear rigidity and Coulomb
    dynamic friction law.

  Thermal Non-rotational Elastic
    Linear elastic repulsion as well as heat transfer.

  Thermal Rotational Bonded
    Same as "Rotational Bonded" with addition of heat transfer.

  Thermal Rotational Friction
    Same as "Rotational Friction" with addition of heat transfer and
    heat generation during frictional slip.

Fixed objects
-------------

Particles not only interact with other particles, but also with
*fixed* objects within the model. These fixed objects are not
subject to the laws of motion and provide a means of imposing
particular types of boundary conditions. Fixed objects include:

  Walls
    An infinite plane characterized by position and normal direction.
  Linear Mesh
    A piecewise linear mesh which can be used to represent a surface in 2D.
  Triangular mesh
    A triangular mesh which can be used to represent a surface in 3D.


"""

__citSection = \
"""
References
==========

.. [CundallGeotech1979] P.A. Cundall and O.A.D Strack
   (1979)
   "A Discrete Numerical Model for Granular Assemblies",
   *Ge\'otechnique*,
   **vol. 29**,
   pp. 47-65.

.. [MoraPAGEOPH1994] P. Mora and D. Place
   (1994)
   "Simulation of the Stick-Slip Instability",
   *Pure Appl. Geophys.*,
   **vol. 143**,
   pp. 61-87.

.. [PlaceJCompPhys1999] D. Place and P. Mora
   (1999)
   "The Lattice Solid Model to Simulate the Physics of Rocks and Earthquakes:
   Incorporation of Friction",
   *J. Comp. Physics*,
   **vol. 150**,
   pp. 332-372.



"""

__doc__ = \
      esys.lsm.doc.Util.setSectionDoc("ESySParticleOverviewSection",__esysParticleOverviewSection) \
    + esys.lsm.doc.Util.setSectionDoc("CitationSection",__citSection) \
    + ("\n:summary: {0:s} overview.\n".format(esys.lsm.util.InstallInfo.pkgName))
