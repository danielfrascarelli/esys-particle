#!/bin/sh
pdflatex ESyS-ParticleNotes
bibtex ESyS-ParticleNotes
./makeglos.pl ESyS-ParticleNotes.glo
pdflatex ESyS-ParticleNotes
pdflatex ESyS-ParticleNotes
#evince ESyS-ParticleNotes.pdf
mv ESyS-ParticleNotes.pdf ../
rm *.aux *.bbl *.blg *.log *.out *.toc
