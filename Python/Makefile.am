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

SUBDIRS = MpiPython BoostPythonUtil esys epydoc .

##
## Hook for building the python API HTML documentation
## Define it using the install-exec-hook because all the binaries
## and python source/byte-code files need to be installed
## so that epydoc can import modules and access the docstrings.
##
if EPYDOC_ENABLED

##
## Need to prefix paths with DESTDIR so that staging works
## correctly (ie make DESTDIR=/tmp/stage install)
##
LSMDATADOCDIR   = $(DESTDIR)$(lsmdatadocdir)
PYTHONAPIDOCDIR = $(DESTDIR)$(pythonapidocdir)
ESYSDIR         = $(DESTDIR)$(esysdir)
LSMDIR          = $(DESTDIR)$(lsmdir)
LSMVISDIR       = $(DESTDIR)$(lsmvisdir)
LIBDIR          = $(DESTDIR)$(libdir)

epydoc: install-exec-hook

install-exec-hook:
	export LD_LIBRARY_PATH=$(LIBDIR):$$LD_LIBRARY_PATH; \
	export PYTHONPATH=$(ESYSDIR):$$PYTHONPATH; \
	if test -d "$(LSMDATADOCDIR)"; then chmod 755 "$(LSMDATADOCDIR)"; else ${INSTALL} -m 755 -d "$(LSMDATADOCDIR)"; fi; \
	if test -d "$(PYTHONAPIDOCDIR)"; then chmod 755 "$(PYTHONAPIDOCDIR)"; else ${INSTALL} -m 755 -d "$(PYTHONAPIDOCDIR)"; fi; \
	if test -f "$(LSMDIR)/__init__.py"; then cd $(PYTHONAPIDOCDIR) && \
		$(EPYDOC_EXE) --verbose --html --no-sourcecode --top "esys.lsm.LsmPy" --css accessBlue.css --name "$(PACKAGE_NAME) $(VERSION)" --url "$(PACKAGE_HOMEPAGE)" `${FIND} $(ESYSDIR) -name *.py` `${FIND} $(ESYSDIR) -name [^lib]*Py$(PYTHON_DLL_EXT)`; fi

uninstall-local:
	if test -d "$(PYTHONAPIDOCDIR)"; then chmod -R +rw "$(PYTHONAPIDOCDIR)" &&\
		rm -r -f "$(PYTHONAPIDOCDIR)" ; fi

endif
