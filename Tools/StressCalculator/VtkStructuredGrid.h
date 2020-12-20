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


#ifndef ESYS_LSM_VTKSTRUCTUREDGRID_H
#define ESYS_LSM_VTKSTRUCTUREDGRID_H

#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "Foundation/vec3.h"
#include "Geometry/Vec3L.h"

#include "Tools/StressCalculator/VtkPiece.h"

namespace esys
{
  namespace lsm
  {
    namespace vtk
    {
      template <typename TmplPointType, typename TmplPointDataTypeTuple>
      class StructuredPiece : public Piece<TmplPointType, TmplPointDataTypeTuple>
      {
      public:
        typedef Piece<TmplPointType, TmplPointDataTypeTuple> Inherited;
        typedef typename Inherited::PointType          PointType;
        typedef typename Inherited::PointValue         PointValue;        
        typedef typename Inherited::PointDataTypeTuple PointDataTypeTuple;
        typedef typename Inherited::PointData          PointData;

        StructuredPiece(const PointType &pointType, const PointDataTypeTuple &pointDataType)
          : Inherited(pointType, pointDataType),
            m_minExtent(),
            m_maxExtent()
        {
        }

        virtual ~StructuredPiece()
        {
        }

        void setExtent(const Vec3L &minIndex, const Vec3L &maxIndex)
        {
          m_minExtent = minIndex;
          m_maxExtent = maxIndex;
        }

        virtual void writeXml(std::ostream &oStream)
        {
          oStream
            << "<Piece Extent=\""
            << getMinExtent()[0] << " "
            << getMaxExtent()[0] << " "
            << getMinExtent()[1] << " "
            << getMaxExtent()[1] << " "
            << getMinExtent()[2] << " "
            << getMaxExtent()[2] << "\">"
            << std::endl;
          this->writePointsXml(oStream);
          this->writePointDataXml(oStream);
          this->writeCellDataXml(oStream);
          oStream << "</Piece>";
        }

        const Vec3L &getMinExtent() const
        {
          return m_minExtent;
        }

        const Vec3L &getMaxExtent() const
        {
          return m_maxExtent;
        }

      protected:

      private:
        Vec3L m_minExtent;
        Vec3L m_maxExtent;
      };

      class StructuredGrid
      {
      private:
        typedef std::vector<XmlPiece *> PiecePtrVector;

      public:
        StructuredGrid()
          : m_pieceVector(),
            m_minExtent(),
            m_maxExtent()
        {
        }

        virtual ~StructuredGrid()
        {
        }

        void setExtent(const Vec3L &minIndex, const Vec3L &maxIndex)
        {
          m_minExtent = minIndex;
          m_maxExtent = maxIndex;
        }

        template <typename TmplStructuredPiece>
        void addPiece(TmplStructuredPiece &piece)
        {
          XmlPiece &xmlPiece = dynamic_cast<XmlPiece &>(piece);
          m_pieceVector.push_back(&xmlPiece);
        }

        virtual void writeXml(std::ostream &oStream)
        {
          oStream 
            << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n"
            << "<StructuredGrid WholeExtent=\""
            << m_minExtent.X() << " "
            << m_maxExtent.X() << " "
            << m_minExtent.Y() << " "
            << m_maxExtent.Y() << " "
            << m_minExtent.Z() << " "
            << m_maxExtent.Z() << "\">"            
            << std::endl;
          for (
            PiecePtrVector::const_iterator it = m_pieceVector.begin();
            it != m_pieceVector.end();
            it++
          )
          {
            (*it)->writeXml(oStream);
            oStream << "\n";
          }
          oStream << "</StructuredGrid>\n";
          oStream << "</VTKFile>";
        }
      private:
        PiecePtrVector m_pieceVector;      
        Vec3L           m_minExtent;
        Vec3L           m_maxExtent;
      };
    }
  }
}

#endif
