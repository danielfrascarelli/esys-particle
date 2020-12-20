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

import vtk as kwvtk

class EdgeExtractor(core.EdgeExtractor):
    def __init__(
        self,
        pointListMap    = lambda dataRecord: dataRecord.getPointList(),
        radiusListMap   = lambda dataRecord: 1.0,
        colorValListMap = lambda dataRecord: 0,
        radiusScale     = 1.0
    ):
        core.EdgeExtractor.__init__(
            self,
            pointListMap,
            radiusListMap,
            colorValListMap,
            radiusScale
        )

    def getResizedList(self, lst, newLen):
        if (not hasattr(lst, "__iter__")):
            newLst = [lst]
        else:
            newLst = list(lst)
        if (len(newLst) < newLen):
            newLst = newLst*((newLen//(len(newLst))) + 1)
        del newLst[newLen:]
        return newLst

    def getVtkIdList(self, idList):
        vtkIdList = kwvtk.vtkIdList()
        for id in idList:
            vtkIdList.InsertNextId(id)
        return vtkIdList

    def getVtkData(self, data):
        pts = kwvtk.vtkPoints()
        cls = kwvtk.vtkCellArray()
        radiusArray = kwvtk.vtkDoubleArray()
        radiusArray.SetName(self.getRadiusFieldName())
        radiusArray.SetNumberOfComponents(1)
        colorArray  = kwvtk.vtkDoubleArray()
        colorArray.SetName(self.getColorFieldName())
        colorArray.SetNumberOfComponents(1)

        ptIndex = 0
        for record in data:
            ptList = self.getPointList(record)
            numPts = len(ptList)
            radiusList = \
              self.getResizedList(self.getRadiusList(record), numPts)
            colorList  = \
              self.getResizedList(self.getColorValList(record), numPts)
            for i in range(0, numPts):
                pts.InsertNextPoint(ptList[i])
                radiusArray.InsertNextValue(radiusList[i]*self.getRadiusScale())
                colorArray.InsertNextValue(colorList[i])
            cls.InsertNextCell(self.getVtkIdList(list(range(ptIndex,ptIndex+numPts))))
            ptIndex += numPts

        vtkPolyData = kwvtk.vtkPolyData()
        vtkPolyData.SetPoints(pts)
        vtkPolyData.SetLines(cls)
        vtkPolyData.GetPointData().SetScalars(radiusArray)
        vtkPolyData.GetPointData().AddArray(colorArray)
        vtkPolyData.GetPointData().SetActiveScalars(radiusArray.GetName())

        return vtkPolyData

    def getColorFieldName(self):
        return "pointDataColorField"

    def getColorFieldComponentIndex(self):
        return 0

    def getRadiusFieldName(self):
        return "pointDataRadiusField"

    def getRadiusFieldComponentIndex(self):
        return 0
