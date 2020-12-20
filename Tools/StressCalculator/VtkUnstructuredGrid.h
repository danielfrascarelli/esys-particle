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


#ifndef ESYS_LSM_VTKUNSTRUCTUREDGRID_H
#define ESYS_LSM_VTKUNSTRUCTUREDGRID_H

#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "Tools/StressCalculator/VtkDataType.h"
#include "Tools/StressCalculator/VtkDataArray.h"
#include "Tools/StressCalculator/VtkDataTypeTuple.h"
#include "Tools/StressCalculator/VtkPiece.h"

namespace esys
{
  namespace lsm
  {
    namespace vtk
    {
      template <typename TmplPointType, typename TmplPointDataTypeTuple>
      class UnstructuredPiece : public Piece<TmplPointType, TmplPointDataTypeTuple>
      {
      public:
        typedef Piece<TmplPointType, TmplPointDataTypeTuple> Inherited;
        typedef typename Inherited::PointType          PointType;
        typedef typename Inherited::PointValue         PointValue;        
        typedef typename Inherited::PointDataTypeTuple PointDataTypeTuple;
        typedef typename Inherited::PointData          PointData;

        UnstructuredPiece(const PointType &pointType, const PointDataTypeTuple &pointDataType)
          : Inherited(pointType, pointDataType)
        {
        }

        virtual void writeXml(std::ostream &oStream)
        {
          oStream
            << "<Piece NumberOfPoints=" << quote(this->getNumPoints())
            << " NumberOfCells=" << quote(this->getNumCells()) << ">" << std::endl;
          this->writePointsXml(oStream);
          this->writePointDataXml(oStream);
          this->writeCellsXml(oStream);
          this->writeCellDataXml(oStream);

          oStream << "</Piece>";
        }

      private:
      };

      class UnstructuredGrid
      {
      private:
        typedef std::vector<XmlPiece *> PiecePtrVector;
      
      public:
        UnstructuredGrid()
          : m_pieceVector()
        {
        }

        virtual ~UnstructuredGrid()
        {
        }

        void addPiece(XmlPiece &piece)
        {
          m_pieceVector.push_back(&piece);
        }
        
        virtual void writeXml(std::ostream &oStream)
        {
          oStream 
            << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n"
            << "<UnstructuredGrid>" << std::endl;
          for (
            PiecePtrVector::const_iterator it = m_pieceVector.begin();
            it != m_pieceVector.end();
            it++
          )
          {
            (*it)->writeXml(oStream);
            oStream << "\n";
          }
          oStream << "</UnstructuredGrid>\n";
          oStream << "</VTKFile>";
        }
      private:
        PiecePtrVector m_pieceVector;
      };
    }
  }
}

#endif
