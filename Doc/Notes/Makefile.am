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

NOTESPDF = ESyS-Particle_Notes.pdf
NOTESNAME = ESyS-ParticleNotes
NOTESSOURCE = \
	ESyS-ParticleNotes.tex \
	ReportStyle.sty \

EXTRA_DIST = $(NOTESSOURCE)

docdir = @lsmdatadocdir@

if PDFDOC_ENABLED

MAKEGLOS = $(srcdir)/makeglos.pl
NOTESAUX = $(NOTESNAME).aux bodytext.aux examples/*.aux tables/*.aux

CLEANFILES = $(NOTESNAME).pdf $(NOTESNAME).log $(NOTESNAME).idx $(NOTESNAME).out \
             $(NOTESNAME).toc $(NOTESAUX) $(NOTESNAME).bbl $(NOTESNAME).blg \
			 $(NOTESNAME).lof $(NOTESPDF)

$(NOTESPDF): $(NOTESSOURCE)
	test "x$(srcdir)" != "x." && cp -pr $(srcdir)/* .
	test -d examples || mkdir examples
	test -d tables || mkdir tables
	$(PDFLATEX) $(NOTESNAME)
	$(BIBTEX) $(NOTESNAME) || true
	$(MAKEGLOS) $(NOTESNAME).glo
	$(PDFLATEX) $(NOTESNAME)
	$(PDFLATEX) $(NOTESNAME)
	cp $(NOTESNAME).pdf $(NOTESPDF)

doc_DATA = $(NOTESPDF)

else

doc_DATA = $(srcdir)/../$(NOTESPDF)

endif