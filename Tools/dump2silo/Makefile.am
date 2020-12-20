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

include $(top_srcdir)/Config/Cpp.am

bin_PROGRAMS = dump2silo

AM_CPPFLAGS += $(SILO_CPPFLAGS) -I$(top_srcdir)/Foundation -I$(top_srcdir)/Geometry

AM_LDFLAGS += $(SILO_LDFLAGS)
 
dump2silo_LDADD = \
    $(top_builddir)/Foundation/libFoundation.la \
    $(top_builddir)/Geometry/libGgGeometry.la \
    $(top_builddir)/Parallel/libParallel.la

dump2silo_SOURCES = \
	main.cpp \
	frame_silo.h \
	frame_silo.cpp 
 
