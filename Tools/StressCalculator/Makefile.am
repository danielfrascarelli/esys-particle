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

lib_LTLIBRARIES = libStressCalculator.la

AM_CPPFLAGS += -I$(top_srcdir)

libStressCalculator_la_LIBADD =                          \
    $(top_builddir)/Foundation/libFoundation.la

libStressCalculator_la_SOURCES = Raw2InteractionReader.h Raw2InteractionReader.cpp InteractionToStressConverter.h InteractionToStressConverter.cpp ContactCollection.cpp ContactCollection.h Contact.h StressTensorCalculator.cpp StressTensorCalculator.h StressTensorCollection.cpp StressTensorCollection.h StressTensor.h VtkUnstructuredGrid.cpp VtkUnstructuredGrid.h VtkDataTypeTuple.h VtkDataTypeTuple.cpp VtkDataArray.h VtkDataType.h CartesianGrid.cpp CartesianGrid.h Vec3Comparer.h VtkStructuredGrid.cpp VtkStructuredGrid.h GaussianGridder.cpp GaussianGridder.h VtkPiece.cpp VtkPiece.h EigenvalueCalculator.cpp EigenvalueCalculator.h Contact.cpp

bin_PROGRAMS = raw2tostress
raw2tostress_SOURCES =                 \
    Main.cpp

raw2tostress_LDADD =                   \
    $(BOOST_FILESYSTEM_LIBS) \
    libStressCalculator.la \
    $(top_builddir)/Foundation/libFoundation.la

