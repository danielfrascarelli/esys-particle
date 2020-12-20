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

bin_PROGRAMS = dump2pov

AM_CPPFLAGS += -I$(top_srcdir)/Foundation -I$(top_srcdir)/Geometry


dump2pov_LDADD = \
    $(top_builddir)/Foundation/libFoundation.la \
    $(top_builddir)/Geometry/libGgGeometry.la \
    $(top_builddir)/Model/libModel.la \
    $(top_builddir)/Parallel/libParallel.la \
    $(top_builddir)/tml/message/libTmlMessage.la \
    $(top_builddir)/tml/type/libTmlType.la

dump2pov_SOURCES = \
	colormap.cpp  \
	colormap3.cpp \
	frame.cpp \
	mesh.cpp \
	geocolormap.cpp\
	main.cpp \
	campos.cpp \
	colormap.h \
	colormap3.h \
	frame.h \
	geocolormap.h \
	mesh.h \
	campos.h 
