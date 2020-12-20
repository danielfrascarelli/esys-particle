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

PYTHON_LTLDFLAGS = -L$(PYTHON_SYS_PREFIX)/lib -lpython$(PYTHON_VERSION)$(PYTHON_VERSION_SUFFIX)

bin_SCRIPTS =             \
    mpipython

EXTRA_DIST = $(bin_SCRIPTS)

bin_PROGRAMS=esysparticle

esysparticle_LDADD=                     \
    $(top_builddir)/Foundation/libFoundation.la \
    $(top_builddir)/Parallel/libParallel.la \
    $(top_builddir)/Geometry/libGgGeometry.la \
    $(top_builddir)/Model/libModel.la \
    $(top_builddir)/tml/message/libTmlMessage.la \
    $(top_builddir)/tml/type/libTmlType.la \
    $(BOOST_FILESYSTEM_LIBS)

esysparticle_LDFLAGS = $(AM_LDFLAGS) $(PYTHON_LTLDFLAGS)

esysparticle_SOURCES=MpiPythonMain.cpp

AM_CPPFLAGS += \
    -I$(top_builddir) \
    -I$(top_srcdir)/Foundation \
    -I$(top_srcdir)/Model      
