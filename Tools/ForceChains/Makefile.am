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

bin_PROGRAMS = fcconv

AM_CPPFLAGS += -I$(top_srcdir)/Foundation


fcconv_LDADD = \
    $(top_builddir)/Foundation/libFoundation.la


fcconv_SOURCES = \
	main.cpp \
	vvf.cpp  \
	vvf.h    \
	vtk.cpp  \
	vtk.h
 


