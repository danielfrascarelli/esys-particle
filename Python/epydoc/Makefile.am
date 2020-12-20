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

SUBDIRS = .

EPYDOCSTYLESHEET    = accessBlue.css
EPYDOCSRCSTYLESHEET = $(srcdir)/$(EPYDOCSTYLESHEET)

GIFFILES = $(srcdir)/AuScopeBannerCropped739x70.gif

EXTRA_DIST = $(GIFFILES) $(EPYDOCSTYLESHEET)
##
## Hook for installing images and stylesheet for the
## the python epydoc API HTML documentation
##
if EPYDOC_ENABLED

##
## Need to prefix paths with DESTDIR so that staging works
## correctly (ie make DESTDIR=/tmp/stage install)
##
LSMDATADOCDIR    = $(DESTDIR)$(lsmdatadocdir)
PYTHONAPIDOCDIR  = $(DESTDIR)$(pythonapidocdir)
EPYDOCTOPDIR     = $(PYTHONAPIDOCDIR)/html
EPYDOCPUBLICDIR  = $(EPYDOCTOPDIR)/public
EPYDOCPRIVATEDIR = $(EPYDOCTOPDIR)/private

install-exec-hook:
	if test -d "$(EPYDOCTOPDIR)"; then chmod 755 "$(EPYDOCTOPDIR)"; else ${INSTALL} -m 755 -d "$(EPYDOCTOPDIR)"; fi; \
	if test -d "$(EPYDOCPUBLICDIR)"; then chmod 755 "$(EPYDOCPUBLICDIR)"; else ${INSTALL} -m 755 -d "$(EPYDOCPUBLICDIR)"; fi; \
	if test -d "$(EPYDOCPRIVATEDIR)"; then chmod 755 "$(EPYDOCPRIVATEDIR)"; else ${INSTALL} -m 755 -d "$(EPYDOCPRIVATEDIR)"; fi; \
	${INSTALL} -m 644 $(EPYDOCSRCSTYLESHEET) $(PYTHONAPIDOCDIR) ; \
	${INSTALL} -m 644 $(GIFFILES) $(EPYDOCTOPDIR) ; \
	${INSTALL} -m 644 $(GIFFILES) $(EPYDOCPUBLICDIR) ; \
	${INSTALL} -m 644 $(GIFFILES) $(EPYDOCPRIVATEDIR)

endif

