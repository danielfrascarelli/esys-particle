#!/bin/sh
pdflatex paper
bibtex paper
./makeglos.pl paper.glo
pdflatex paper
pdflatex paper
#evince paper.pdf
cp paper.pdf ../ESyS-Particle_Tutorial.pdf
