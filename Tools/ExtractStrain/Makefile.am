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

bin_PROGRAMS = strainextract

AM_CPPFLAGS += -I$(top_srcdir)/Foundation -I$(top_srcdir)/ntable/src


strainextract_LDADD = \
    $(top_builddir)/Foundation/libFoundation.la \
    $(top_builddir)/ntable/src/libntable.la


strainextract_SOURCES = \
	main.cpp \
	DataParticle.cpp \
	DataExtractor.cpp \
	SnapFileHelp.cpp
	 


