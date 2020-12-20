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

lib_LTLIBRARIES = libBoostPythonUtil.la

libBoostPythonUtil_la_SOURCES = ListConverter.h ListConverter.cpp Util.h \
	PythonIterIterator.h PythonIterIterator.hpp
libBoostPythonUtil_la_LIBADD = $(BOOST_PYTHON_LIBS)

