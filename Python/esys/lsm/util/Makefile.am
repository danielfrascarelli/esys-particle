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

if VTK_ENABLED
  HAVE_VTK_BOOLEAN=True
else
  HAVE_VTK_BOOLEAN=False
endif

if POVRAY_ENABLED
  HAVE_POVRAY_BOOLEAN=True
else
  HAVE_POVRAY_BOOLEAN=False
endif

lsmutil_LTLIBRARIES = libFoundationPy.la FoundationPy.la

libFoundationPy_la_LDFLAGS = -shared
libFoundationPy_la_SOURCES = Vec3Py.cpp \
	Vec3Py.h \
	BoundingBoxPy.cpp \
	BoundingBoxPy.h \
	BoundingSpherePy.cpp \
	BoundingSpherePy.h \
	RngPy.cpp \
	RngPy.h \
	SetPy.cpp \
	SetPy.h \
	SetPy.hpp \
	VectorPy.cpp \
	VectorPy.h \
	VectorPy.hpp \
	QuaternionPy.cpp \
	QuaternionPy.h 	

libFoundationPy_la_LIBADD  =                                     \
    $(top_builddir)/Foundation/libFoundation.la                  \
    $(top_builddir)/Python/BoostPythonUtil/libBoostPythonUtil.la \
    $(BOOST_PYTHON_LIBS)

FoundationPy_la_LDFLAGS = -module -shared -shrext $(PYTHON_DLL_EXT)
FoundationPy_la_SOURCES = ExportModuleUtil.cpp
FoundationPy_la_LIBADD  = libFoundationPy.la


dist_lsmutil_PYTHON =   \
    __init__.py         \
    FileUtil.py         \
    OptParse.py         \
    pathSearcher.py     \
    process.py

nodist_lsmutil_PYTHON = \
    InstallInfo.py

CLEANFILES = InstallInfo.py

EXTRA_DIST=InstallInfo.tmpl

InstallInfo.py: InstallInfo.tmpl Makefile
	$(PYTHON) -c 'from __future__ import print_function;\
	  import os;\
	  writefile=open("InstallInfo.py", "w");\
	  readfile=open(os.path.join("$(srcdir)", "InstallInfo.tmpl"), "r");\
	  text=str.join("",(readfile.readlines())).format\
	    (\
	      "$(prefix)",\
	      "$(pythondir)",\
	      "$(PACKAGE_NAME)",\
	      "$(VERSION)",\
	      "$(PACKAGE_HOMEPAGE)",\
	      "$(pkgdatadir)",\
	      "$(POVRAY_EXE)",\
	      "$(HAVE_VTK_BOOLEAN)",\
	      "$(HAVE_POVRAY_BOOLEAN)"\
	    );\
	  print(text,file=writefile);\
	  readfile.close();\
	  writefile.close();\
	  '

