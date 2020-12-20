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

from esys.lsm.vis import core

import types
import vtk as kwvtk
import sys

"""
From http://www.python.org/dev/peps/pep-0404/: "Python 2 has two basic
integer types: a native machine-sized "int" type, and an arbitrary-length
"long" type. These have been merged in Python 3 into a single "int" type
analogous to Python 2's "long" type."

This test preserves the behavior of code written for the "long" type under
Python 2.
"""
if sys.version_info[0] > 2:
  lint = int
else:
  lint = long
  
typeVtkDataArrayDict = dict()
typeVtkDataArrayDict[int]   = kwvtk.vtkIntArray
typeVtkDataArrayDict[float] = kwvtk.vtkDoubleArray
typeVtkDataArrayDict[lint]  = kwvtk.vtkLongArray

def getVtkDataArrayClass(thing):
    if (hasattr(thing, "__len__")):
        return typeVtkDataArrayDict[type(thing[0])]
    return typeVtkDataArrayDict[type(thing)]

class FieldMap:
    def __init__(self, name, map):
        self.name = name
        self.map  = map
        self.vtkDataArray = None
        self.componentIdxList = []

    def initialise(self, record):
        val = self.map(record)
        self.vtkDataArray = getVtkDataArrayClass(val)()
        self.vtkDataArray.SetName(self.name)
        self.vtkDataArray.SetNumberOfComponents(getattr(val, "__len__", lambda:1)())
        self.componentIdxList = list(range(0,self.vtkDataArray.GetNumberOfComponents()))

    def insertNext(self, record):
        tupleIdx = self.vtkDataArray.GetNumberOfTuples()
        tupyl = self.map(record)
        for j in self.componentIdxList:
            self.vtkDataArray.InsertComponent(tupleIdx, j, tupyl[j])

    def getVtkDataArray(self):
        return self.vtkDataArray

class GlyphExtractor:
    def __init__(self, pointsMap, fieldMapList):
        self.pointsMap = pointsMap
        self.fieldMapList = fieldMapList

    def addRecordToVtkData(self, record, vtkPoints):
        vtkPoints.InsertNextPoint(self.pointsMap(record))
        for fieldMap in self.fieldMapList:
            fieldMap.insertNext(record)

    def getVtkData(self, data):
        rIt = iter(data)
        record = next(rIt)
        for fieldMap in self.fieldMapList:
            fieldMap.initialise(record)

        vtkPoints = kwvtk.vtkPoints()
        self.addRecordToVtkData(record, vtkPoints)
        for record in rIt:
            self.addRecordToVtkData(record, vtkPoints)

        vtkUnstrGrid = kwvtk.vtkUnstructuredGrid()
        vtkUnstrGrid.SetPoints(vtkPoints)
        for fieldMap in self.fieldMapList:
            vtkUnstrGrid.GetPointData().AddArray(fieldMap.getVtkDataArray())
        return vtkUnstrGrid

