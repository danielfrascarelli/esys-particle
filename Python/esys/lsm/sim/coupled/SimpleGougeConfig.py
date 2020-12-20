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

from __future__   import division
from esys.escript import *
from esys.finley import Rectangle,Brick
from esys.finley import Merge

from esys.lsm.geometry import *
from esys.lsm.util     import *

from .MeshToVtkConverter import MeshToVtkConverter

def roundToInt(flt):
    return int(round(flt))

class SimpleGougeMesh:
    def __init__(
        self,
        boundingBox1,
        boundingBox2,
        gridSize,
    ):
        self.bBox1         = boundingBox1
        self.bBox2         = boundingBox2
        self.gridSize      = gridSize
        self.mesh          = None

    def getDim(self):
        dim = 3
        if (self.bBox1.getMinPt()[2] == self.bBox1.getMaxPt()[2]):
           dim = 2
        return dim

    def getGridPointsPerDim(self):
        sizes = (self.bBox1.getSize(), self.bBox2.getSize())
        return \
            (
               list(map(roundToInt, (sizes[0]/self.gridSize).toList())),
               list(map(roundToInt, (sizes[1]/self.gridSize).toList()))
            )

    def generate(self):
        (numGridPts1, numGridPts2) = self.getGridPointsPerDim()
        periodic0 = True
        if (self.getDim() == 2):
            mesh1 = \
                Rectangle(
                    numGridPts1[0],
                    numGridPts1[1],
                    order=1,
                    l0=self.bBox1.getSize()[0],
                    l1=self.bBox1.getSize()[1],
                    periodic0=periodic0
                )
            mesh2 = \
                Rectangle(
                    numGridPts2[0],
                    numGridPts2[1],
                    order=1,
                    l0=self.bBox2.getSize()[0],
                    l1=self.bBox2.getSize()[1],
                    periodic0=periodic0
                )
        else:
            mesh1 = \
                Brick(
                    numGridPts1[0],
                    numGridPts1[1],
                    numGridPts1[2],
                    order=1,
                    l0=self.bBox1.getSize()[0],
                    l1=self.bBox1.getSize()[1],
                    l2=self.bBox1.getSize()[2],
                    periodic0=periodic0
                )
            mesh2 = \
                Brick(
                    numGridPts2[0],
                    numGridPts2[1],
                    numGridPts2[2],
                    order=1,
                    l0=self.bBox2.getSize()[0],
                    l1=self.bBox2.getSize()[1],
                    l2=self.bBox2.getSize()[2],
                    periodic0=periodic0
                )
        x1 = mesh1.getX()
        x2 = mesh2.getX()
        #
        # To get the same identical tags on the particle boundary elements
        # of both meshes, reflect about the y-axis.
        # 
        x2[1] = -x2[1]
        x2[1] += abs(inf(x2[1]))

        #
        # Now translate the meshes so they correspond to the
        # specified bounding boxes.
        #
        for d in range(0, self.getDim()):
            x1[d] += self.bBox1.getMinPt()[d]
            x2[d] += self.bBox2.getMinPt()[d]

        mesh1.setX(x1)
        mesh2.setX(x2)
        self.mesh = Merge([mesh1, mesh2])

    def write(self, fileName):
        self.mesh.write(fileName)

class SimpleGougeConfig:
    def __init__(
        self,
        particleGougePrms,
        continuumBlockHeight,
        fileNamePrefix="gouge",
        gridSize = None
    ):
        self.gougeBlock = GougeBlock(particleGougePrms)
        if (gridSize == None):
            gridSize = 2.0*particleGougePrms.getMinRadius()
        self.gridSize = gridSize
        self.continuumHeight = continuumBlockHeight
        self.fileNamePrefix = fileNamePrefix
        self.meshGenerator = None

    def generateDomain(self):
        self.gougeBlock.generate()
        particleBBox = self.gougeBlock.getDomainBoundingBox()
        minPt = particleBBox.getMinPt()
        maxPt = particleBBox.getMaxPt()
        upperMinPt = Vec3(minPt[0], maxPt[1], minPt[2])
        upperMaxPt = maxPt + Vec3(0, self.continuumHeight, 0)
        lowerMinPt = minPt - Vec3(0, self.continuumHeight, 0)
        lowerMaxPt = Vec3(maxPt[0], minPt[1], maxPt[2])

        self.meshGenerator = \
            SimpleGougeMesh(
                BoundingBox(lowerMinPt, lowerMaxPt),
                BoundingBox(upperMinPt, upperMaxPt),
                self.gridSize
            )
        self.meshGenerator.generate()

    def getParticleGeometryXmlFileName(self):
        return self.fileNamePrefix + "ParticlesVtk.xml"

    def getMeshXmlFileName(self):
        return self.fileNamePrefix + "MeshVtk.xml"

    def getParticleGeometryFileName(self):
        return self.fileNamePrefix + ".geo"

    def getMeshFileName(self):
        return self.fileNamePrefix + ".msh"

    def writeVtkXml(self):
        self.write()
        self.gougeBlock.writeVtkXml(self.getParticleGeometryXmlFileName())
        MeshToVtkConverter(self.getMeshFileName(), self.getMeshXmlFileName()).convert()

    def write(self):
        self.gougeBlock.write(self.getParticleGeometryFileName())
        self.meshGenerator.write(self.getMeshFileName())

