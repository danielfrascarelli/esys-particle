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
# Init file for esys.lsm package.

from . import Logging
from .util  import InstallInfo
from esys.lsm.LsmPy import *

__doc__ = \
"""
Lattice Solid Model discrete element simulation package. For an introduction, 
see the tutorial at the U{{{pkgName} home page <{pkgHomePageUrl}>}}.
@organization: The University of Queensland Centre for Geoscience Computing\
 (Brisbane, Queensland, Australia) U{{earth.uq.edu.au/centre-geoscience-computing}}
@copyright: (C) 2003-2017 by The University of Queensland Centre for Geoscience Computing
@license: Open Software License version 3.0\
 U{{www.apache.org/licenses/LICENSE-2.0}}
@undocumented: Logging
""".format(pkgHomePageUrl=InstallInfo.pkgHomePageUrl, pkgName=InstallInfo.pkgName)

LsmPy.__doc__ = \
"""
Dynamic linked library module defining classes for creating
particle models - this library exposes C++ implementations in Python.
@group Model: LsmMpi,*Sphere*,Runnable
@group Interaction Types: *InteractionGroup,*InteractionPrms,*GravityPrms, *BuoyancyPrms,*DampingPrms,*ElasticPrms,*FrictionPrms,*BondPrms, *MeshPrms,*WallPrms,*BeamPrms,*DashpotPrms,*TriggerPrms,*2DPrms
@group Saving Data: *SaverPrms
@group Saving Model State: *CheckPointPrms
@group Miscellaneous: __*,set*,check*,*BuildPrms,ParticleIdPair*
"""

checkMpiDimensions.__doc__ = \
"""
Checks validity of C{numProcesses} and C{mpiDimList} arguments.
C{numProcesses} should be positive, C{mpiDimList} should contain
3 elements and the product of non-zero elements should be
a divisor of C{numProcesses}. Raises C{ValueError} if arguments
are determined to be invalid.
@type numProcesses: int
@kwarg numProcesses: number of mpi worker processes.
@type mpiDimList: list of three ints
@kwarg mpiDimList: MPI cartesian communicator list.
@raises ValueError: if C{numProcesses} and C{mpiDimList} are not consistent.
"""

checkParticleType.__doc__ = \
"""
Raises C{ValueError} if particle type is invalid.
@type particleType: string
@kwarg particleType: Type of particle.
"""
