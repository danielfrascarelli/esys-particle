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

lsmsimcoupled_LTLIBRARIES = libLsmEscriptPy.la LsmEscriptPy.la

libLsmEscriptPy_la_LDFLAGS = -shared
libLsmEscriptPy_la_SOURCES =                                     \
     LsmMpiEscriptPy.cpp LsmMpiEscriptPy.h

libLsmEscriptPy_la_LIBADD  =                                     \
    $(top_builddir)/Python/BoostPythonUtil/libBoostPythonUtil.la

LsmEscriptPy_la_LDFLAGS = -module -shared -shrext $(PYTHON_DLL_EXT)
LsmEscriptPy_la_SOURCES = ExportModuleLsmEscript.cpp
LsmEscriptPy_la_LIBADD  =                        \
    $(top_builddir)/Python/esys/lsm/libLsmPy.la  \
    libLsmEscriptPy.la

lsmsimcoupled_PYTHON =    \
    __init__.py           \
    SimpleGouge.py        \
    SimpleGougeConfig.py  \
    MeshToVtkConverter.py \
    DxToVtkConverter.py

