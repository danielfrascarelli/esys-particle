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

bin_PROGRAMS = fracextract

AM_CPPFLAGS += -I$(top_srcdir)/Foundation -I$(top_srcdir)/Geometry


fracextract_LDADD = \
    $(top_builddir)/Foundation/libFoundation.la \
    $(top_builddir)/Parallel/libParallel.la \
    $(top_builddir)/Geometry/libGgGeometry.la

fracextract_SOURCES = \
	fracframe.cpp  \
	fracframe.h \
	fracwriter.cpp \
	fracwriter.h \
	FitPlane.h \
	FitPlane.cpp \
	frac_dist.h \
	frac_dist.cpp \
	main.cpp   

