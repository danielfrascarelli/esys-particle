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

import os
import re

VTK_LINE  = 3
VTK_PIXEL = 8
VTK_VOXEL = 11

class ParseError(Exception):
    pass

def findFileMatch(regExprString, f):
    expr = re.compile(regExprString)
    line = f.readline()
    match = expr.match(line)
    while ((match == None) and (len(line) > 0)):
        line = f.readline()
        match = expr.match(line)
    return match

class Node:
    def __init__(self, mshLine=None, ptLine=None):
        self.id = None
        self.idCirc = None
        self.tag = None
        self.pt  = [0.0]*3
        
        if (mshLine != None):
            self.initFromMeshLine(mshLine)
        elif (ptLine != None):
            self.initFromPtLine(ptLine)

    def initFromMeshLine(self, line):
        elemList = str.split(str.strip(line), " ")
        self.id = int(elemList[0])
        self.idCirc = int(elemList[1])
        self.tag = int(elemList[2])
        ptList  = list(map(float, elemList[3:]))
        self.pt  = [0.0]*3
        self.pt[0:len(ptList)] = ptList

    def initFromPtLine(self, line):
        elemList = str.split(str.strip(line), " ")
        ptList  = list(map(float, elemList))
        self.pt  = [0.0]*3
        self.pt[0:len(ptList)] = ptList

    def translateBy(self, displacement):
        self.pt[0:len(displacement)] = \
          [self.pt[i] + displacement[i] for i in range(0, len(displacement))]

    def writePt(self, f):
        f.write(str.join(" ", (list(map(str, self.pt)))))

numNodesVtkCellTypeDict = dict()
numNodesVtkCellTypeDict[2] = VTK_LINE
numNodesVtkCellTypeDict[4] = VTK_PIXEL
numNodesVtkCellTypeDict[8] = VTK_VOXEL

class Cell:
    def __init__(self, meshLine=None, nodeIdLine=None):
        self.id     = None
        self.tag    = None
        self.nodeIdList = []
        
        if (meshLine != None):
            self.initFromMeshLine(meshLine)
        elif (nodeIdLine != None):
            self.initFromNodeIdLine(nodeIdLine)

    def initFromMeshLine(self, line):
        elemList    = str.split(str.strip(line), " ")
        self.id     = int(elemList[0])
        self.tag    = int(elemList[1])
        self.nodeIdList = list(map(int, elemList[2:]))

    def initFromNodeIdLine(self, line):
        elemList    = str.split(str.strip(line), " ")
        self.nodeIdList = list(map(int, elemList))
    
    def getIndexMapList(self):
        indexMapList = list(range(0, len(self.nodeIdList)))
        if (len(self.nodeIdList) == 4):
            indexMapList = [0, 1, 3, 2]
        if (len(self.nodeIdList) == 8):
            indexMapList = [0, 1, 3, 2, 4, 5, 7, 6]
        return indexMapList

    def getNumNodes(self):
        return len(self.nodeIdList)

    def getVtkCellType(self):
        return numNodesVtkCellTypeDict[self.getNumNodes()]

    def writeNodeIdList(self, f):
        nodeIdList = [self.nodeIdList[i] for i in self.getIndexMapList()]
        f.write(str.join("", (map(str, nodeIdList))))

class NodeList:
    def __init__(self):
        self.f = None
        self.nodeList = []

    def __len__(self):
        return len(self.nodeList)

    def __getitem__(self, i):
        return self.nodeList[i]

    def __setitem__(self, i, node):
        self.nodeList[i] = node

    def findMatch(self, regExprString):
        return findFileMatch(regExprString, self.f)

    def getDim(self):
        return self.dim

    def findNodes(self):
        match = self.findMatch("^(\d)D-Nodes\s*(\d*)")
        if (match != None):
            self.dim =      int(match.group(1))
            self.numNodes = int(match.group(2))
        else:
            raise ParseError("No nodes header found in file " + self.fileName)

    def append(self, node):
        self.nodeList.append(node)

    def readNodes(self):
        self.findNodes()
        for i in range(0, self.numNodes):
            self.append(Node(self.f.readline()))

    def read(self, f):
        self.f = f
        self.readNodes()
        self.f = None

    def writeVtkPointsDataArray(self, f):
        f.write("<DataArray NumberOfComponents=\"{0:d}\" type=\"Float32\" format=\"ascii\">\n".format(3))
        for node in self.nodeList:
            node.writePt(f)
            f.write("\n")
        f.write("</DataArray>\n")

class CellList:
    def __init__(self):
        self.f = None
        self.cellType     = None
        self.nodesPerCell = None
        self.cellList = []

    def append(self, cell):
        self.cellList.append(cell)

    def __len__(self):
        return len(self.cellList)

    def __getitem__(self, i):
        return self.cellList[i]

    def __setitem__(self, i, cell):
        self.cellList[i] = cell

    def findCells(self):
        match = findFileMatch("^(\D+)(\d+)\s*(\d+)", self.f)
        if (match != None):
            self.cellType     = match.group(1)
            self.nodesPerCell = int(match.group(2))
            self.numCells     = int(match.group(3))
        else:
            raise ParseError("No cells header found in file " + self.fileName)

    def append(self, cell):
        self.cellList.append(cell)

    def readCells(self):
        self.findCells()
        for i in range(0, self.numCells):
            self.append(Cell(self.f.readline()))

    def read(self, f):
        self.f = f
        self.readCells()
        self.f = None

    def writeVtkCellsDataArray(self, f):
        self.writeVtkConnectivity(f)
        self.writeVtkOffsets(f)
        self.writeVtkTypes(f)

    def writeVtkConnectivity(self, f):
        f.write(
            "<DataArray NumberOfComponents=\"1\"" +\
            " type=\"Int32\" Name=\"connectivity\"" +\
            " format=\"ascii\">\n"
        )
        for cell in self.cellList:
            cell.writeNodeIdList(f)
            f.write("\n")
        f.write("</DataArray>\n")

    def writeVtkOffsets(self, f):
        f.write(
            "<DataArray NumberOfComponents=\"1\"" +\
            " type=\"Int32\" Name=\"offsets\"" +\
            " format=\"ascii\">\n"
        )
        offset = 0
        for cell in self.cellList:
            offset += cell.getNumNodes()
            f.write(str(offset))
            f.write("\n")
        f.write("</DataArray>\n")

    def getVtkCellType(self):
        if (self.cellType == "Line"):
            return VTK_LINE
        elif (self.cellType == "Rec"):
            return VTK_PIXEL
        elif (self.cellType == "Hex"):
            return VTK_VOXEL

    def writeVtkTypes(self, f):
        f.write(
            "<DataArray NumberOfComponents=\"1\"" +\
            " type=\"Int32\" Name=\"types\"" +\
            " format=\"ascii\">\n"
        )
        for cell in self.cellList:
            f.write(str(cell.getVtkCellType()))
            f.write("\n")
        f.write("</DataArray>\n")
 
class UnstructuredGrid:
    def __init__(self, nodeList, cellList):
        self.nodeList = nodeList
        self.cellList = cellList

    def getNumNodes(self):
        return len(self.nodeList)

    def getNumCells(self):
        return len(self.cellList)

    def writeVtkXml(self, fileName):
        f = file(fileName, "w")
        f.write("<?xml version=\"1.0\"?>\n")
        f.write("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n")
        f.write("<UnstructuredGrid>\n")
        f.write(
            "<Piece NumberOfPoints=\"{0:d}\" NumberOfCells=\"{1:d}\">\n"\
            .format(self.getNumNodes(), self.getNumCells())
        )
        f.write("<Points>\n")
        self.nodeList.writeVtkPointsDataArray(f)
        f.write("</Points>\n")
        f.write("<PointData>\n")
        f.write("</PointData>\n")
        f.write("<Cells>\n")
        self.cellList.writeVtkCellsDataArray(f)
        f.write("</Cells>\n")
        f.write("<CellData>\n")
        f.write("</CellData>\n")
        f.write("</Piece>\n")
        f.write("</UnstructuredGrid>\n")
        f.write("</VTKFile>\n")

class MeshReader:
    def __init__(self, fileName):
        self.fileName = fileName
        self.f = None
        self.dim = None
        self.nodeList = NodeList()
        self.cellList = CellList()
        self.bndryCellList = CellList()

    def getNumNodes(self):
        return len(self.nodeList)

    def getNumCells(self):
        return len(self.cellList)

    def getNumBndryCells(self):
        return len(self.bndryCellList)

    def read(self):
        self.f = file(self.fileName, "r")
        self.nodeList.read(self.f)
        self.cellList.read(self.f)
        self.bndryCellList.read(self.f)
        self.f.close()
        self.f = None

    def writeVtkXml(self, fileName):
        UnstructuredGrid(self.nodeList, self.cellList).writeVtkXml(fileName)

class MeshToVtkConverter:
    def __init__(self, meshFileName, vtkXmlFileName):
        self.meshFileName = meshFileName
        self.vtkXmlFileName = vtkXmlFileName

    def convert(self):
        meshData = MeshReader(self.meshFileName)
        meshData.read()
        meshData.writeVtkXml(self.vtkXmlFileName)
