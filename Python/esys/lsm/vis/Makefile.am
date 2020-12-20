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

include $(top_srcdir)/Config/Python.am

if VTK_ENABLED
  VTK     = vtk
endif

if POVRAY_ENABLED
  POVRAY     = povray
endif

SUBDIRS = core $(VTK) $(POVRAY) 

lsmvis_PYTHON = __init__.py

AM_CPPFLAGS += $(PYTHON_CPPFLAGS)

