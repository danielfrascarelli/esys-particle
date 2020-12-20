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

lib_LTLIBRARIES=libTmlComm.la

libTmlComm_la_SOURCES = cart_comm.cpp comm.cpp comm_world.cpp cart_comm.h cart_comm.hpp comm_coll.hpp comm.h comm.hpp comm_world.h

libTmlComm_la_LIBADD =                                    \
	$(top_builddir)/Foundation/libFoundation.la       \
	$(top_builddir)/tml/message/libTmlMessage.la      \
	$(top_builddir)/tml/type/libTmlType.la
