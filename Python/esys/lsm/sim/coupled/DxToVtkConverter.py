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

from . import MeshToVtkConverter
from .MeshToVtkConverter import CellList
from .MeshToVtkConverter import Node, NodeList
from .MeshToVtkConverter import UnstructuredGrid, ParseError

import re

class Cell(MeshToVtkConverter.Cell):
    def getIndexMapList(self):
        return list(range(0, len(self.nodeIdList)))

class DxObject:
    def __init__(self, id, dxClass, type, rank, shape, items):
        self.id = id
        self.dxClass = dxClass
        self.type = type
        self.rank = rank
        self.shape = shape
        self.items = items
    
    def getNumItems(self):
        return self.items

class DxData:
    def __init__(self):
        self.cellList = CellList()
        self.nodeList = NodeList()
        self.f        = None
    
    def findNextObjectLine(self):
        regex = re.compile(
            "\s*object\s+(\d+)\s+class\s+(.+)\s+type\s+(.+)\s+rank\s+(\d+)\s+" +\
            "shape\s+(\d+)\s+items\s+(\d+)\s+data\s+follows"
        )
        line = self.f.readline()
        match = regex.match(line)
        while ((len(line) > 0) and (match == None)):
            line = self.f.readline()
            match = regex.match(line)
        if (match != None):
            return \
                DxObject(
                    id = int(match.group(1)),
                    dxClass = match.group(2),
                    type = match.group(3),
                    rank = int(match.group(4)),
                    shape = int(match.group(5)),
                    items = int(match.group(6))
                )
        else:
            raise ParseError("Object match not found.")

    def readNodes(self):
        dxObject = self.findNextObjectLine()
        for i in range(0, dxObject.getNumItems()):
            self.nodeList.append(Node(ptLine=self.f.readline()))

    def readCells(self):
        dxObject = self.findNextObjectLine()
        for i in range(0, dxObject.getNumItems()):
            self.cellList.append(Cell(nodeIdLine=self.f.readline()))

    def getPosition(self, line):
        return list(map(float, str.split(str.strip(line), " ")))

    def translateNodes(self):
        dxObject = self.findNextObjectLine()
        for i in range(0, dxObject.getNumItems()):
            self.nodeList[i].translateBy(self.getPosition(self.f.readline()))

    def read(self, fileName):
        self.f = file(fileName, "r")
        self.readNodes()
        self.readCells()
        self.translateNodes()

    def writeVtkXml(self, fileName):
        UnstructuredGrid(self.nodeList, self.cellList).writeVtkXml(fileName)
