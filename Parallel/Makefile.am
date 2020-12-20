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

lib_LTLIBRARIES=libParallel.la

AM_CPPFLAGS += \
      $(SILO_CPPFLAGS) \
    -I$(top_builddir) \
    -I$(top_srcdir)/Foundation \
    -I$(top_srcdir)/Model      \
    -I$(top_srcdir)/Fields     \
    -I$(top_srcdir)/Interface  \
    -I$(top_srcdir)/MainBE     \
    -I$(top_srcdir)/Socket

AM_LDFLAGS += $(SILO_LDFLAGS)

libParallel_la_LDFLAGS =

libParallel_la_SOURCES = \
    ASubLattice.cpp \
		ASubLattice.h \
		LatticeMaster.cpp \
		LatticeMaster.h \
		mpia2abuf.cpp \
		mpia2abuf.h \
		mpibarrier.cpp \
		mpibarrier.h \
		mpibuf.cpp \
		mpibuf.h \
		mpicmdbuf.cpp \
		mpicmdbuf.h \
		mpisgbuf.cpp \
		mpisgbuf.h \
		mpisgvbuf.cpp \
		mpisgvbuf.h \
		mpi_tag_defs.h \
		mpivbuf.cpp \
		mpivbuf.h \
		sublattice_cmd.h \
		SubLatticeControler.cpp \
		SubLatticeControler.h \
		SubLattice.cpp \
		SubLattice.h \
		SubLattice.hpp \
		RotSubLattice.h \
		RotSubLattice.hpp \
		CheckPointer.cpp \
	 	CheckPointController.cpp \
	 	CheckPointParams.cpp \
		MpiInfo.cpp \
		CheckPointer.h \
		CheckPointController.h \
	 	CheckPointParams.h \
		MpiInfo.h \
		GeometryReader.cpp \
		GeometryReader.h \
		CheckPointInfo.cpp \
		CheckPointLoader.cpp \
		CheckPointInfo.h \
	 	CheckPointLoader.h \
		LatticeParam.cpp \
		LatticeParam.h \
	 	MpiWrap.h \
		BroadCast_cmd.h \
		BroadCast_cmd.hpp \
	 	BroadCast_cmd.cpp \
	 	BMesh2D_cmd.h \
		BMesh2D_cmd.cpp \
	 	GetRef_cmd.h \
		GetRef_cmd.cpp \
		RankAndComm.h \
		LatticeMaster.hpp \
		TagData.h \
		TagData.cpp \ 	
		GMRESSolverMaster.h \
		GMRESSolverMaster.cpp \
		GMRESSolverSlave.h \
		GMRESSolverSlave.cpp

libParallel_la_LIBADD =                                \
    $(top_builddir)/Foundation/libFoundation.la        \
    $(top_builddir)/Geometry/libGgGeometry.la          \
    $(top_builddir)/ppa/src/libppa.la                  \
    $(top_builddir)/tml/message/libTmlMessage.la       \
    $(top_builddir)/tml/type/libTmlType.la             \
    $(top_builddir)/Model/libModel.la                  \
    $(top_builddir)/tml/comm/libTmlComm.la             \
    $(top_builddir)/Fields/libFields.la \
    $(top_builddir)/pis/libpis.la

noinst_HEADERS = \
	MeshReaderImpl.h \
	Mesh2DReaderImpl.h  
