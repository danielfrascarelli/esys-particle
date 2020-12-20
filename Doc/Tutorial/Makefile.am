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

TUTPDF = ESyS-Particle_Tutorial.pdf
TUTNAME = paper
TUTSOURCE = \
	paper.tex \
	bodytext.tex \
	epsf.sty \
	glossary.sty \
	bingle_output.jpg \
	bingle_vis_snaps_2.jpg \
	bingle_vis_snaps.jpg \
	gravity_cube_snaps.jpg \
	gravity_plot.jpg \
	hopper_flow.jpg \
	PV_examples.png \
	PV_initial_window.png \
	rotational_bonds.jpg \
	SC_bulk_friction.png \
	SC_dilatancy.png \
	ShearCellSetup.jpg \
	slope_fail_snaps.jpg \
	slope_friction_quarter_snaps.jpg \
	slope_friction_snaps.jpg \
	slope_friction_wall_snaps.jpg \
	stress_strain_diagram.jpg \
	trimesh_normal.jpg \
	UCS_epot.jpg \
	UCS_force.jpg \
	UCS_nbonds.jpg \
	UCS_stress_strain.jpg \
	examples/bingle_chk.py.tex \
	examples/bingle_output.py.tex \
	examples/bingle.py.tex \
	examples/bingle_Runnable.py.tex \
	examples/bingle_vis.py.tex \
	examples/floorMesh.msh.tex \
	examples/gravity_cube.py.tex \
	examples/gravity.py.tex \
	examples/hopper_flow.py.tex \
	examples/POVsnaps.py.tex \
	examples/rot_compress.py.tex \
	examples/ServoWallLoader.py.tex \
	examples/shearcell.py.tex \
	examples/slope_fail.py.tex \
	examples/slope_friction_floor.py.tex \
	examples/slope_friction.py.tex \
	examples/slope_friction_walls.py.tex \
	examples/WallLoader.py.tex \
	tables/field_names.tex \
	tables/file_formats.tex \
	tables/interaction_groups.tex


EXTRA_DIST = $(TUTSOURCE)

docdir = @lsmdatadocdir@

if PDFDOC_ENABLED

MAKEGLOS = $(srcdir)/makeglos.pl
TUTAUX = $(TUTNAME).aux bodytext.aux examples/*.aux tables/*.aux

CLEANFILES = $(TUTNAME).pdf $(TUTNAME).log $(TUTNAME).idx $(TUTNAME).out \
             $(TUTNAME).toc $(TUTAUX) $(TUTNAME).bbl $(TUTNAME).blg \
			 $(TUTNAME).lof $(TUTPDF)

$(TUTPDF): $(TUTSOURCE)
	##
	## If building within the source directory, this test fails,
	## and so make exits and the build also fails.  To handle a 
	## possible test failure from 'test "x$(srcdir)" != "x."',
	## frame it within 'if... then... fi...' as in the next
	## statement.
	##
	##test "x$(srcdir)" != "x." && cp -pr $(srcdir)/* .
	if test "x$(srcdir)" != "x."; then cp -pr $(srcdir)/* .; fi;
	test -d examples || mkdir examples
	test -d tables || mkdir tables
	$(PDFLATEX) $(TUTNAME)
	$(BIBTEX) $(TUTNAME) || true
	$(MAKEGLOS) $(TUTNAME).glo
	$(PDFLATEX) $(TUTNAME)
	$(PDFLATEX) $(TUTNAME)
	cp $(TUTNAME).pdf $(TUTPDF)

doc_DATA = $(TUTPDF)

else

doc_DATA = $(srcdir)/../$(TUTPDF)

endif

