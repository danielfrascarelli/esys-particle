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


#ifndef ESYS_LSM_BPULISTCONVERTER_H
#define ESYS_LSM_BPULISTCONVERTER_H

#include <boost/python.hpp>

#include "Python/BoostPythonUtil/Util.h"

#include <vector>

namespace esys
{
  namespace lsm
  {
    namespace bpu
    {
      template <typename TmplValueType>
      class DefaultExtractor
      {
      public:
        typedef TmplValueType value_type;
        
        value_type operator()(const boost::python::object &pyObject) const
        {
          return boost::python::extract<value_type>(pyObject);
        }
      };

      template<typename TmplValue>
      std::vector<TmplValue> listToVector(const boost::python::list &pythonList)
      {
        DefaultExtractor<TmplValue> extractor;
        std::vector<TmplValue> vec;
        const int numElements = esys::lsm::bpu::len(pythonList);
        vec.reserve(numElements);
        for (int i = 0; i < numElements; i++)
        {
          vec.push_back(extractor(pythonList[i]));
        }
        return vec;
      }

      template<typename TmplValue>
      std::vector<TmplValue> tupleToVector(const boost::python::tuple &pythonTulple)
      {
        DefaultExtractor<TmplValue> extractor;
        std::vector<TmplValue> vec;
        const int numElements = esys::lsm::bpu::len(pythonTulple);
        vec.reserve(numElements);
        for (int i = 0; i < numElements; i++)
        {
          vec.push_back(extractor(pythonTulple[i]));
        }
        return vec;
      }

      template<typename TmplValue, typename TmplExtractor>
      std::vector<TmplValue> listToVector(const boost::python::list &pythonList, TmplExtractor extractor=TmplExtractor())
      {
        std::vector<TmplValue> vec;
        const int numElements = esys::lsm::bpu::len(pythonList);
        vec.reserve(numElements);
        for (int i = 0; i < numElements; i++)
        {
          vec.push_back(extractor(pythonList[i]));
        }
        return vec;
      }

      template <typename TmplVector>
      boost::python::list vectorToList(const TmplVector &vec)
      {
        boost::python::list pythonList;
        for (
          typename TmplVector::const_iterator it = vec.begin();
          it != vec.end();
          it++
        )
        {
          pythonList.append(*it);
        }
        return pythonList;
      }
    }
  }
}

#endif
