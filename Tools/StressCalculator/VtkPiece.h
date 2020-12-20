/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2017 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0              //
//                                                         //
/////////////////////////////////////////////////////////////


#ifndef ESYS_LSM_VTKPIECE_H
#define ESYS_LSM_VTKPIECE_H

#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "Tools/StressCalculator/VtkDataType.h"
#include "Tools/StressCalculator/VtkDataArray.h"
#include "Tools/StressCalculator/VtkDataTypeTuple.h"

namespace esys
{
  namespace lsm
  {
    namespace vtk
    {
      class XmlPiece
      {
      public:
        virtual void writeXml(std::ostream &oStream) = 0;
      };

      template <typename TmplPointType, typename TmplPointDataTypeTuple>
      class Piece : public XmlPiece
      {
      public:
        typedef TmplPointType                               PointType;
        typedef typename PointType::value_type              PointValue;        
        typedef TmplPointDataTypeTuple                      PointDataTypeTuple;
        typedef typename PointDataTypeTuple::DataValueTuple PointData;

        Piece(const PointType &pointType, const PointDataTypeTuple &pointDataType)
          : m_pointData(pointDataType),
            m_pointValueArray(pointType),
            m_pointIndexMap()
        {
        }

        virtual ~Piece()
        {
        }

        int getIndex(const PointValue &point) const
        {
          typename PointIndexMap::const_iterator it = m_pointIndexMap.find(point);
          return ((it != m_pointIndexMap.end()) ? it->second : -1);
        }

        void setPoint(const PointValue &point, const PointData &data)
        {
          int index = getIndex(point);
          if (index < 0) {
            index = m_pointValueArray.size();
            m_pointIndexMap.insert(
              typename PointIndexMap::value_type(point, index)
            );
            m_pointValueArray.setData(index, point);
          }
          m_pointData.setData(index, data);
        }
        
        int getNumPoints() const
        {
          return m_pointValueArray.size();
        }

        int getNumCells() const
        {
          return 0;
        }
        
        virtual void writeXml(std::ostream &oStream) = 0;
        
        virtual void writePointsXml(std::ostream &oStream)
        {
          oStream << "<Points>" << "\n";
          m_pointValueArray.writeXml(oStream);
          oStream << "\n</Points>" << "\n";
        }

        virtual void writePointDataXml(std::ostream &oStream)
        {
          oStream << "<PointData>" << "\n";
          m_pointData.writeXml(oStream);
          oStream << "\n</PointData>" << "\n";
        }

        virtual void writeCellsXml(std::ostream &oStream)
        {
          oStream << "<Cells>" << "\n";
          DataArray<Int32Type> connectivity(Int32Type("connectivity"));
          connectivity.writeXml(oStream);
          oStream << "\n";
          DataArray<Int32Type> offsets(Int32Type("offsets"));
          offsets.writeXml(oStream);
          oStream << "\n";
          DataArray<UInt8Type> types(UInt8Type("types"));
          types.writeXml(oStream);
          oStream << "\n</Cells>" << "\n";
        }

        virtual void writeCellDataXml(std::ostream &oStream)
        {
          oStream << "<CellData>" << "\n";
          oStream << "</CellData>" << "\n";        
        }

      private:
        typedef DataArray<PointType>      PointValueArray;
        typedef std::map<PointValue, int> PointIndexMap;

        PointDataTypeTuple m_pointData;
        PointValueArray    m_pointValueArray;
        PointIndexMap      m_pointIndexMap;
      };
    }
  }
}

#endif
