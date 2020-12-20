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

##
include $(top_srcdir)/Config/Python.am

SUBDIRS =        \
    util         \
    .            \
    geometry     \
    vis          

lsm_PYTHON = __init__.py Logging.py 

AM_CPPFLAGS += $(SILO_CPPFLAGS)
AM_LDFLAGS += $(SILO_LDFLAGS)

lsm_LTLIBRARIES = libLsmPy.la LsmPy.la

libLsmPy_la_LDFLAGS = -shared
libLsmPy_la_SOURCES = LsmMpiPy.cpp \
	LsmMpiPy.h \
	RunnablePy.cpp \
	RunnablePy.h \
	CheckPointParamsPy.cpp \
	CheckPointParamsPy.h \
	InteractionParamsPy.cpp \
	InteractionParamsPy.h \
	InteractionGroupPy.cpp \
	InteractionGroupPy.h \
	BondInteractionGroupPy.cpp \
	BondInteractionGroupPy.h \
	BondedTriMeshPrmsPy.cpp \
	BondedTriMeshPrmsPy.h \
	BondedMesh2DPrmsPy.cpp \
	BondedMesh2DPrmsPy.h \
	MeshBuildParamsPy.cpp \
	MeshBuildParamsPy.h \
	WallPrmsPy.cpp \
	WallPrmsPy.h \
	SphereBodyPrmsPy.cpp \
	SphereBodyPrmsPy.h \
	ParticlePy.cpp \
	ParticlePy.h \
	FieldSaverPrmsPy.cpp \
	FieldSaverPrmsPy.h \
	ParticleFieldSaverPrmsPy.cpp \
	ParticleFieldSaverPrmsPy.h \
	InteractionFieldSaverPrmsPy.cpp \
	InteractionFieldSaverPrmsPy.h \
	ElasticTriMeshPrmsPy.cpp \
	ElasticTriMeshPrmsPy.h \
	WallFieldSaverPrmsPy.cpp \
	WallFieldSaverPrmsPy.h \
	TriangleFieldSaverPrmsPy.cpp \
	TriangleFieldSaverPrmsPy.h \
	RotParticlePy.cpp \
	RotParticlePy.h \
	LmParticleAdder.cpp \
	LmParticleAdder.h \
	LmParticleAdder.hpp \
	TriggerPrmsPy.cpp \
	TriggerPrmsPy.h \
	ElasticMesh2DPrmsPy.cpp \
	ElasticMesh2DPrmsPy.h \
	ParticleIdPairSetPy.cpp \
	ParticleIdPairSetPy.h \
	ParticleIdPairPy.cpp \
	ParticleIdPairPy.h \
	ParticleIdPairVectorPy.cpp \
	ParticleIdPairVectorPy.h \
	RotParticleViPy.cpp \
	RotParticleViPy.h \
	RotThermalParticlePy.cpp \
	RotThermalParticlePy.h \
	FluidFieldSaverPrmsPy.cpp \
	FluidFieldSaverPrmsPy.h \
	FluidInteractionFieldSaverPrmsPy.cpp \
	FluidInteractionFieldSaverPrmsPy.h

libLsmPy_la_LIBADD  =                                            \
    $(top_builddir)/Foundation/libFoundation.la                  \
    $(top_builddir)/Model/libModel.la                            \
    $(top_builddir)/Parallel/libParallel.la                      \
    $(top_builddir)/Python/BoostPythonUtil/libBoostPythonUtil.la \
    $(top_builddir)/Python/esys/lsm/util/libFoundationPy.la

LsmPy_la_LDFLAGS = -module -shared -shrext $(PYTHON_DLL_EXT)
LsmPy_la_SOURCES = ExportModuleLsm.cpp
LsmPy_la_LIBADD  = libLsmPy.la                

