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

include $(top_srcdir)/Config/Cpp.am

#SUBDIRS = dump2geo dump2pov mesh2pov rotextract pfcconvert StressCalculator MacroFrictionCalculator ExtractGrains dump2vtk ExtractProfile ExtractFractures ForceChains slip2vtk ExtractStrain
SUBDIRS = dump2geo dump2pov mesh2pov rotextract StressCalculator ExtractGrains dump2vtk ExtractFractures ForceChains ExtractStrain MergeLog


if SILO_ENABLED
SUBDIRS += dump2silo
endif

