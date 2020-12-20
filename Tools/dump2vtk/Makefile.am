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

## Makefile.am for dump2vtk
include $(top_srcdir)/Config/Cpp.am

bin_PROGRAMS = dump2vtk

AM_CPPFLAGS += -I$(top_srcdir)/Foundation -I$(top_srcdir)/Geometry


dump2vtk_LDADD = \
    $(top_builddir)/Foundation/libFoundation.la \
    $(top_builddir)/Geometry/libGgGeometry.la \
    $(top_builddir)/Parallel/libParallel.la \
    $(BOOST_FILESYSTEM_LIBS)

dump2vtk_SOURCES = \
	main.cpp \
	frame_vtk.h \
	frame_vtk.cpp 
 
