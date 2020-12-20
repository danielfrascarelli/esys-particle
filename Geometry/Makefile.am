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

lib_LTLIBRARIES = libGgGeometry.la

AM_CPPFLAGS += -I$(top_srcdir)/Foundation

libGgGeometry_la_SOURCES = ARandomAssembly.cpp ARandomAssembly.h \
					BasicInteraction.cpp BasicInteraction.h LayeredBlock.cpp LayeredBlock.h Line.cpp Line.h \
					PaddedBlock3D.cpp PaddedBlock3D.h Plane3D.cpp Plane3D.h RandomAssembly2D.cpp \
					RandomAssembly2D.h RandomAssembly3D.cpp RandomAssembly3D.h RandomBlock3D.cpp \
					RandomBlock3D.h RandomBlock.cpp RandomBlock.h SimpleNTable3D.cpp SimpleNTable3D.h \
					SimpleNTable.cpp SimpleNTable.h Sphere2d.cpp Sphere2d.h Sphere3d.cpp Sphere3d.h Sphere.h \
					SplitBlock3D.cpp SplitBlock3D.h SplitBlock.cpp SplitBlock.h GougeBlock3D.cpp \
					ParticleGenerator.cpp RegularBlockGenerator.cpp RandomBlockGenerator.cpp BlockGenerator.cpp \
					GeometryInfo.cpp GeometryInfo.h Vec3L.cpp VtkXmlWriter.cpp ParticleFitter.cpp \
					GranularGougeBlock3D.cpp SphAggGougeBlock.cpp SphAggGougeBlock.h SphAggGougeBlock.hpp \
					SphereBlockGenerator.cpp Triangle.cpp Triangle.h AEdge.h AEdge.cpp Edge.h Edge.cpp Edge2D.h \
					Edge2D.cpp Corner.h Corner.cpp Corner2D.h Corner2D.cpp SphereBoxVolCalculator.cpp \
					SphereBoxVolCalculator.h CircleBoxVolCalculator.cpp CircleBoxVolCalculator.h \
					IntersectionVolCalculator.cpp IntersectionVolCalculator.h IntersectionVolCalculator.hpp \
					SimpleBlock.cpp SimpleBlock.h SimpleBlock.hpp DistConnections.cpp DistConnections.h \
					DistConnections.hpp ParticleCollection.cpp ParticleCollection.h ParticleCollection.hpp \
					CubicBlock.cpp CubicBlock.h CubicBlock.hpp GridIterator.cpp GridIterator.h Grain.cpp \
					Grain.h Grain.hpp GrainCollection.cpp GrainCollection.h GrainCollection.hpp \
					Packer.cpp Packer.h Packer.hpp BoxPacker.cpp BoxPacker.h BoxPacker.hpp \
					CubicBoxPacker.cpp CubicBoxPacker.h CubicBoxPacker.hpp RandomBoxPacker.cpp \
					RandomBoxPacker.h RandomBoxPacker.hpp SphereFitter.cpp SphereFitter.h SphereFitter.hpp \
					CircularNeighbourTable.cpp CircularNeighbourTable.h NeighbourTable.cpp NeighbourTable.h \
					GrainRandomBoxPacker.cpp GrainRandomBoxPacker.h GrainRandomBoxPacker.hpp SimpleParticleData.cpp \
					SimpleParticleData.h SimpleParticleData.hpp SimpleParticle.cpp SimpleParticle.h \
					SimpleParticle.hpp PackerGenerators.cpp PackerGenerators.h PackerGenerators.hpp \
			RandomSpherePacker.cpp RandomSpherePacker.h RandomSpherePacker.hpp ParticleComparer.cpp \
			ParticleComparer.h ParticleComparer.hpp GougeConfig.cpp GougeConfig.h GougeConfig.hpp \
			PointCloud.h 	PointCloud.cpp CubicBlockIterator.cpp CubicBlockIterator.h \
		CubicBlockIterator.hpp ClosePackIterator.cpp ClosePackIterator.h ClosePackIterator.hpp \
		HexagBlockIterator.cpp HexagBlockIterator.h HexagBlockIterator.hpp ClosePackBlock.cpp \
		ClosePackBlock.h ClosePackBlock.hpp HexagBlock.cpp HexagBlock.h HexagBlock.hpp \
		ClosePackOrientation.h CircularNeighbourTable.hpp NeighbourTable.hpp \
		LineSegment.cpp LineSegment.h FaultedBlock2d.cpp FaultedBlock2d.h RoughPaddedBlock3d.h RoughPaddedBlock3d.cpp \
		SphereNeighbours.cpp SphereNeighbours.hpp SphereNeighbours.cpp \
	RectPatch.h RectPatch.cpp

libGgGeometry_la_LIBADD =                          \
	$(top_builddir)/Foundation/libFoundation.la

noinst_HEADERS = \
	GridIterator.h \
	GougeBlock3D.h \
	ParticleGenerator.h \
	RegularBlockGenerator.h \
	RandomBlockGenerator.h \
	BlockGenerator.h \
	RandomBlockGenerator.h \
	SimpleParticle.h \
	NeighbourTable.h \
	Vec3L.h \
	CircularNeighbourTable.h \
	VtkXmlWriter.h \
	ParticleFitter.h \
	GranularGougeBlock3D.h \
	SphAggGougeBlock.h \
	SphereBlockGenerator.h
