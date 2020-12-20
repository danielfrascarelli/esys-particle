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

## Makefile.am for dump2geo
include $(top_srcdir)/Config/Cpp.am

bin_PROGRAMS = dump2geo

AM_CPPFLAGS += -I$(top_srcdir)/Foundation -I$(top_srcdir)/Geometry


dump2geo_LDADD = \
    $(top_builddir)/Foundation/libFoundation.la \
    $(top_builddir)/Geometry/libGgGeometry.la \
    $(top_builddir)/Parallel/libParallel.la \
    $(top_builddir)/Model/libModel.la \
    $(top_builddir)/Fields/libFields.la \
    $(top_builddir)/tml/message/libTmlMessage.la

dump2geo_SOURCES = \
	main.cpp \
	frame_geo.h \
	frame_geo.cpp 
 
