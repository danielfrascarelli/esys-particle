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

lib_LTLIBRARIES = libMacroFrictionCalc.la

AM_CPPFLAGS += -I$(top_srcdir)

libMacroFrictionCalc_la_LIBADD =                          \
    $(top_builddir)/Foundation/libFoundation.la

libMacroFrictionCalc_la_SOURCES = MacroFrictionCalculator.cpp MacroFrictionCalculator.h LinearWindowAverager.cpp LinearWindowAverager.h WallForceReader.cpp WallForceReader.h WallForcesToFrictionConverter.cpp WallForcesToFrictionConverter.h

bin_PROGRAMS = forcetofric

forcetofric_SOURCES =                 \
    Main.cpp

forcetofric_LDADD =                   \
    libMacroFrictionCalc.la

