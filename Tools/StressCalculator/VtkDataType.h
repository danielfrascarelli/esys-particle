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


#ifndef ESYS_LSM_VTKDATATYPE_H
#define ESYS_LSM_VTKDATATYPE_H

#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "Foundation/vec3.h"
#include "Foundation/Matrix3.h"

namespace esys
{
  namespace lsm
  {
    namespace vtk
    {
      typedef std::string ValueTypeName;
      typedef std::string FormatTypeName;

      static const ValueTypeName UInt8     = "UInt8";
      static const ValueTypeName Int16     = "Int16";
      static const ValueTypeName UInt16    = "UInt16";
      static const ValueTypeName Int32     = "Int32";
      static const ValueTypeName UInt32    = "UInt32";
      static const ValueTypeName Int64     = "Int64";
      static const ValueTypeName UInt64    = "UInt64";
      static const ValueTypeName Float32   = "Float32";
      static const ValueTypeName Float64   = "Float64";
      
      static const FormatTypeName ascii    = "ascii";
      static const FormatTypeName binary   = "binary";
      static const FormatTypeName appended = "appended";

      template <typename TmplType>
      std::string quote(const TmplType &thing)
      {
        std::stringstream sStream;
        sStream << "\"" << thing << "\"";
        return sStream.str();
      }
      
      template <typename TmplValueType>
      class DataType
      {
      public:
        typedef TmplValueType value_type;
        DataType(
          const ValueTypeName &valueTypeName,
          const std::string &dataName,
          unsigned int numComponents,
          const FormatTypeName &format = ascii,
          unsigned int offset = 0
        ) :
          m_valueTypeName(valueTypeName),
          m_dataName(dataName),
          m_numComponents(numComponents),
          m_format(format),
          m_offset(offset)      
        {
        }
        
        std::string getXmlAttributeString() const
        {
          std::stringstream sStream;
          
          sStream
            << "type=" << quote(m_valueTypeName) << " "
            << "Name=" << quote(m_dataName) << " "
            << "NumberOfComponents=" << quote(m_numComponents) << " "
            << "format=" << quote(m_format);
          if (m_format == appended) {
            sStream << " offset=" << quote(m_offset);
          }
          
          return sStream.str();
        }
        
      private:
        ValueTypeName   m_valueTypeName;
        std::string     m_dataName;
        unsigned int    m_numComponents;
        FormatTypeName  m_format;
        unsigned int    m_offset;
      };
      
      class Float64Type : public DataType<double>
      {
      public:
        typedef DataType<double> Inherited;
        Float64Type(
          const std::string &name,
          const FormatTypeName &format=ascii,
          int offset=0
        ) 
          : Inherited(Float64, name, 1, format, offset)
        {
        }
      };
      
      class Float32Type : public DataType<float>
      {
      public:
        typedef DataType<float> Inherited;
        Float32Type(
          const std::string &name,
          const FormatTypeName &format=ascii,
          int offset=0
        )
          : Inherited(Float32, name, 1, format, offset)
        {
        }
      };

      class UInt8Type : public DataType<unsigned char>
      {
      public:
        typedef DataType<unsigned char> Inherited;
        UInt8Type(
          const std::string &name,
          const FormatTypeName &format=ascii,
          int offset=0
        )
          : Inherited(UInt8, name, 1, format, offset)
        {
        }
      };

      class Int32Type : public DataType<int>
      {
      public:
        typedef DataType<int> Inherited;
        Int32Type(
          const std::string &name,
          const FormatTypeName &format=ascii,
          int offset=0
        )
          : Inherited(Int32, name, 1, format, offset)
        {
        }
      };

      class Vec3Type : public DataType<Vec3>
      {
      public:
        typedef DataType<Vec3> Inherited;
        Vec3Type(
          const std::string &name,
          const FormatTypeName &format=ascii,
          int offset=0
        ) 
          : Inherited(Float64, name, 3, format, offset)
        {
        }
      };
      
      class Matrix3Type : public DataType<Matrix3>
      {
      public:
        typedef DataType<Matrix3> Inherited;
        Matrix3Type(
          const std::string &name,
          const FormatTypeName &format=ascii,
          int offset=0
        ) 
          : Inherited(Float64, name, 9, format, offset)
        {
        }
      };
    }
  }
}

#endif
