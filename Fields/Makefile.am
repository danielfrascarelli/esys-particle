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
include $(top_srcdir)/Config/Cpp.am

lib_LTLIBRARIES=libFields.la

AM_CPPFLAGS += $(SILO_CPPFLAGS) -I$(top_srcdir)/Foundation -I$(top_srcdir)/Model

AM_LDFLAGS += $(SILO_LDFLAGS)
 
libFields_la_SOURCES =\
	CheckedScalarInteractionFieldSlave.h \
	CheckedScalarInteractionFieldSlave.hpp \
	CheckedScalarInteractionFieldSlaveTagged.h \
	CheckedScalarInteractionFieldSlaveTagged.hpp \
	InteractionFieldMaster.cpp \
	InteractionFieldMaster.h \
	InteractionFieldSlave.hpp \
	InteractionFieldSlave.h \
	ParticleFieldMaster.cpp \
	ParticleFieldMaster.h \
	TriggeredVectorParticleFieldMaster.cpp \
	TriggeredVectorParticleFieldMaster.h \
	ScalarParticleDistributionMaster.h \
	ScalarParticleDistributionMaster.cpp \
	ScalarParticleFieldSlaveTagged.h \
	ScalarParticleFieldSlaveTagged.hpp \
	VectorParticleFieldSlaveTagged.h \
	VectorParticleFieldSlaveTagged.hpp \
	ScalarInteractionFieldSlave.h \
	ScalarInteractionFieldSlave.hpp \
	ScalarInteractionFieldSlaveTagged.h \
	ScalarInteractionFieldSlaveTagged.hpp \
	ScalarParticleFieldSlave.h \
	ScalarParticleFieldSlave.hpp \
	ScalarTriangleFieldMaster.h \
	ScalarTriangleFieldMaster.cpp \
	ScalarTriangleFieldSlave.h \
	ScalarTriangleFieldSlave.cpp \
	VectorInteractionFieldSlave.h \
	VectorInteractionFieldSlave.hpp \
	VectorInteractionFieldMaster.h \
	VectorInteractionFieldMaster.cpp \
	VectorParticleFieldSlave.h \
	VectorParticleFieldSlave.hpp \
	VectorTriangleFieldMaster.h \
	VectorTriangleFieldMaster.cpp \
	VectorTriangleFieldSlave.h \
	VectorTriangleFieldSlave.cpp \
	VectorEdge2DFieldSlave.h \
	VectorEdge2DFieldSlave.cpp \
	MaxTrigger.cpp \
	MaxTrigger.h \
	FluidFieldMaster.h \
	FluidFieldMaster.cpp \
	ScalarFluidFieldSlave.hpp \
	ScalarFluidFieldSlave.h \
	VectorFluidFieldSlave.hpp \
	VectorFluidFieldSlave.h \
	FluidInteractionFieldMaster.h \
	FluidInteractionFieldMaster.cpp \
	ScalarFluidInteractionFieldSlave.hpp \
	ScalarFluidInteractionFieldSlave.h \
	VectorFluidInteractionFieldSlave.hpp \
	VectorFluidInteractionFieldSlave.h

libFields_la_LIBADD =                           \
  $(top_builddir)/Foundation/libFoundation.la   \
  $(top_builddir)/tml/type/libTmlType.la        \
  $(top_builddir)/tml/comm/libTmlComm.la        \
  $(top_builddir)/Model/libModel.la             

noinst_HEADERS = 
