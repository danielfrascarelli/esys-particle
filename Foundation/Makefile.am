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

lib_LTLIBRARIES=libFoundation.la

libFoundation_la_LIBADD = $(BOOST_FILESYSTEM_LIBS) $(BOOST_PYTHON_LIBS)
libFoundation_la_LDFLAGS = $(AM_LDFLAGS) $(BOOST_PYTHON_LDFLAGS)
libFoundation_la_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_PYTHON_CPPFLAGS)

libFoundation_la_SOURCES = \
	console.cpp \
	console.hpp \
	console.h \
	Counter.cpp \
	Counter.h \
	Error.h \
	Matrix3.cpp \
	Matrix3.h \
	Matrix3.hpp \
	Quaternion.cpp \
	Quaternion.h \
	Quaternion.hpp \
	t_list.h \
	t_list.hpp \
	vec3.cpp \
	vec3.h \
	vec3.hpp \
	Timer.cpp \
	Timer.h \
	StringUtil.cpp \
	StringUtil.h \
	realdist.h \
	realdist.cpp \
	quadtuple.h \
	StlIterator.h \
	quintuple.h \
	PathSearcher.cpp \
	PathSearcher.h \
	Runnable.cpp \
	Runnable.h \
	PathUtil.cpp \
	PathUtil.h \
	Rng.cpp \
	Rng.h \
	Rng.hpp \
	BoundingSphere.cpp \
	BoundingSphere.h \
	BoundingSphere.hpp \
	BoundingBox.cpp \
	BoundingBox.h \
	BoundingBox.hpp \
	cube_eq.cpp \
	cube_eq.h \
	triplet.h \
	Functional.h

noinst_HEADERS =

