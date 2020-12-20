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

lsmgeometry_LTLIBRARIES = libGeometryPy.la GeometryPy.la

libGeometryPy_la_LDFLAGS = -shared
libGeometryPy_la_SOURCES = SimpleBlockPy.cpp SimpleBlockPy.h MiscPy.cpp \
	MiscPy.h DistConnectionsPy.cpp DistConnectionsPy.h TaggedIdConnectionPy.cpp \
	TaggedIdConnectionPy.h IteratorPy.cpp IteratorPy.h IteratorPy.hpp ParticleCollectionPy.cpp \
	ParticleCollectionPy.h CubicBlockPy.cpp CubicBlockPy.h GrainCollectionPy.cpp \
	GrainCollectionPy.h GrainPy.cpp GrainPy.h PackerPy.cpp PackerPy.h GougeConfigPrmsPy.cpp \
	GougeConfigPrmsPy.h GougeConfigPy.cpp GougeConfigPy.h HexagBlockPy.cpp HexagBlockPy.h \
	SimpleSpherePy.cpp SimpleSpherePy.h \
        SimpleSphereNeighboursPy.cpp SimpleSphereNeighboursPy.h \
	SimpleSphereCollectionPy.h SimpleSphereCollectionPy.cpp
libGeometryPy_la_LIBADD  =                                       \
    $(top_builddir)/Geometry/libGgGeometry.la                    \
    $(top_builddir)/Python/BoostPythonUtil/libBoostPythonUtil.la \
    $(top_builddir)/Python/esys/lsm/libLsmPy.la                  \
    $(top_builddir)/Python/esys/lsm/util/libFoundationPy.la


GeometryPy_la_LDFLAGS = -module -shared -shrext $(PYTHON_DLL_EXT)
GeometryPy_la_SOURCES = ExportModuleGeometry.cpp
GeometryPy_la_LIBADD  = libGeometryPy.la

lsmgeometry_PYTHON =            \
    __init__.py                 \
    SingleParticleGrainGen.py   \
    SphericalGrainGen.py
