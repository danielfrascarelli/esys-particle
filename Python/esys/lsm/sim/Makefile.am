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

AM_CPPFLAGS += $(SILO_CPPFLAGS)
AM_LDFLAGS += $(SILO_LDFLAGS)

lsmsim_PYTHON =        \
    __init__.py        \
    wavePropagation.py

if ESCRIPT_ENABLED
  COUPLED     = coupled
endif

SUBDIRS = $(COUPLED)

##

lsmsim_LTLIBRARIES = libWavePropagationPy.la WavePropagationPy.la

libWavePropagationPy_la_LDFLAGS = -shared
libWavePropagationPy_la_SOURCES =                   \
    WavePropagationPy.cpp WavePropagationPy.h

libWavePropagationPy_la_LIBADD  =                   \
    $(top_builddir)/Python/esys/lsm/libLsmPy.la

WavePropagationPy_la_LDFLAGS = -module -shared -shrext $(PYTHON_DLL_EXT)
WavePropagationPy_la_SOURCES = ExportModuleWavePropagation.cpp
WavePropagationPy_la_LIBADD  = libWavePropagationPy.la

