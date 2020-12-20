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

bin_PROGRAMS = mesh2pov

AM_CPPFLAGS += -I$(top_srcdir)/Foundation -I$(top_srcdir)/Geometry

mesh2pov_LDADD = \
    $(top_builddir)/Foundation/libFoundation.la \
    $(top_builddir)/Geometry/libGgGeometry.la \
    $(top_builddir)/Model/libModel.la \
    $(top_builddir)/Parallel/libParallel.la \
    $(top_builddir)/tml/message/libTmlMessage.la \
    $(top_builddir)/tml/type/libTmlType.la

mesh2pov_SOURCES = \
    main.cpp 

