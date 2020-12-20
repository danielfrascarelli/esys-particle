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


#ifndef ESYS_LSM_VTKDATATYPETUPLE_H
#define ESYS_LSM_VTKDATATYPETUPLE_H

#include <vector>
#include <map>
#include <iostream>

#include <boost/tuple/tuple.hpp>

namespace esys
{
  namespace lsm
  {
    namespace vtk
    {
      template <typename TmplDataType> class DataArray;

      class NullDataType
      {
      public:
        typedef char value_type;
      };

      // Specialisation for the NullDataType DataArray
      template <>
      class DataArray<NullDataType>
      {
      public:
        typedef NullDataType::value_type value_type;
        DataArray(const NullDataType &dataType)
        {
        }

        void setData(int , const value_type &)
        {
        }
        
        void writeXml(std::ostream &)
        {
        }
      };
      
      // Forward declaration
      template <
        typename DT0=NullDataType, typename DT1=NullDataType,
        typename DT2=NullDataType, typename DT3=NullDataType,
        typename DT4=NullDataType, typename DT5=NullDataType,
        typename DT6=NullDataType, typename DT7=NullDataType,
        typename DT8=NullDataType, typename DT9=NullDataType
        >
        class DataTypeTuple;

      // Implementation
      template <typename DT0, typename DT1, typename DT2, typename DT3, typename DT4,
                typename DT5, typename DT6, typename DT7, typename DT8, typename DT9>
      class DataTypeTuple
      {
      public:
        typedef
          boost::tuple<
            typename DT0::value_type, typename DT1::value_type,
            typename DT2::value_type, typename DT3::value_type,
            typename DT4::value_type, typename DT5::value_type,
            typename DT6::value_type, typename DT7::value_type,
            typename DT8::value_type, typename DT9::value_type
          > DataValueTuple;

        typedef DataArray<DT0> DA0;
        typedef DataArray<DT1> DA1;
        typedef DataArray<DT2> DA2;
        typedef DataArray<DT3> DA3;
        typedef DataArray<DT4> DA4;
        typedef DataArray<DT5> DA5;
        typedef DataArray<DT6> DA6;
        typedef DataArray<DT7> DA7;
        typedef DataArray<DT8> DA8;
        typedef DataArray<DT9> DA9;

        typedef
          boost::tuple<
            DA0,DA1,DA2,DA3,DA4,DA5,DA6,DA7,DA8,DA9
          > DataArrayTuple;

        DataTypeTuple(
          const DT0 &dt0 = DT0(),
          const DT1 &dt1 = DT1(),
          const DT2 &dt2 = DT2(),
          const DT3 &dt3 = DT3(),
          const DT4 &dt4 = DT4(),
          const DT5 &dt5 = DT5(),
          const DT6 &dt6 = DT6(),
          const DT7 &dt7 = DT7(),
          const DT8 &dt8 = DT8(),
          const DT9 &dt9 = DT9()
        )
          : m_dataArrayTuple(
              DA0(dt0), DA1(dt1), DA2(dt2), DA3(dt3), DA4(dt4),
              DA5(dt5), DA6(dt6), DA7(dt7), DA8(dt8), DA9(dt9)
            )
        {
        }
        
        void setData(int index, const DataValueTuple &data)
        {
          m_dataArrayTuple.template get<0>().setData(index, data.template get<0>());
          m_dataArrayTuple.template get<1>().setData(index, data.template get<1>());
          m_dataArrayTuple.template get<2>().setData(index, data.template get<2>());
          m_dataArrayTuple.template get<3>().setData(index, data.template get<3>());
          m_dataArrayTuple.template get<4>().setData(index, data.template get<4>());
          m_dataArrayTuple.template get<5>().setData(index, data.template get<5>());
          m_dataArrayTuple.template get<6>().setData(index, data.template get<6>());
          m_dataArrayTuple.template get<7>().setData(index, data.template get<7>());
          m_dataArrayTuple.template get<8>().setData(index, data.template get<8>());
          m_dataArrayTuple.template get<9>().setData(index, data.template get<9>());
        }
        
        void writeXml(std::ostream &oStream)
        {
          m_dataArrayTuple.template get<0>().writeXml(oStream);
          m_dataArrayTuple.template get<1>().writeXml(oStream);
          m_dataArrayTuple.template get<2>().writeXml(oStream);
          m_dataArrayTuple.template get<3>().writeXml(oStream);
          m_dataArrayTuple.template get<4>().writeXml(oStream);
          m_dataArrayTuple.template get<5>().writeXml(oStream);
          m_dataArrayTuple.template get<6>().writeXml(oStream);
          m_dataArrayTuple.template get<7>().writeXml(oStream);
          m_dataArrayTuple.template get<8>().writeXml(oStream);
          m_dataArrayTuple.template get<9>().writeXml(oStream);
        }

      private:
        DataArrayTuple m_dataArrayTuple;
      };
    }
  }
}

#endif
