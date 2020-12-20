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


#ifndef ESYS_LSM_BPUPYTHONITERITERATOR_H
#define ESYS_LSM_BPUPYTHONITERITERATOR_H

#include <boost/python.hpp>
#include <patchlevel.h>

namespace esys
{
  namespace lsm
  {
    namespace bpu
    {
      template <typename TmplExtractType>
      class PythonIterIterator
      {
      public:
        typedef TmplExtractType value_type;
        PythonIterIterator(boost::python::object &iteratable);
  
        bool hasNext() const;
  
        value_type next();
  
        void update();
  
      private:
        bool                  m_hasNext;
        boost::python::object m_next;
        boost::python::object m_iter;
      };
    }
  }
}

#include "Python/BoostPythonUtil/PythonIterIterator.hpp"

#endif
