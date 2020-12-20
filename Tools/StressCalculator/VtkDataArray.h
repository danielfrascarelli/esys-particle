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


#ifndef ESYS_LSM_VTKDATAARRAY_H
#define ESYS_LSM_VTKDATAARRAY_H

#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "Tools/StressCalculator/VtkDataType.h"

namespace esys
{
  namespace lsm
  {
    namespace vtk
    {
      template <typename TmplDataType>
      class DataArray
      {
      public:
        typedef TmplDataType                  DataType;
        typedef typename DataType::value_type value_type;

        DataArray(const DataType &dataType)
          : m_dataType(dataType),
            m_valueVector()
        {
          m_valueVector.reserve(512);
        }
        
        void setData(int index, const value_type &val)
        {
          if (static_cast<int>(m_valueVector.size()) <= index) {
            m_valueVector.resize(index+1);
          }
          m_valueVector.at(index) = val;
        }

        int size() const
        {
          return m_valueVector.size();
        }
        
        std::string getXmlAttributeString() const
        {
          return m_dataType.getXmlAttributeString();
        }
        
        void writeXml(std::ostream &oStream)
        {
          oStream << "<DataArray " << getXmlAttributeString() << ">" << "\n";
          for (typename ValueVector::const_iterator it = m_valueVector.begin(); it != m_valueVector.end(); it++)
          {
            oStream << (*it) << "\n";
          }
          oStream << "</DataArray>";
        }
      private:
        typedef std::vector<value_type> ValueVector;
        DataType    m_dataType;
        ValueVector m_valueVector;
      };
    }
  }
}

#endif
