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


#ifndef ESYS_LSMITERATORPY_H
#define ESYS_LSMITERATORPY_H

#include <boost/python.hpp>
#include <patchlevel.h>
#include <string>

namespace esys
{
  namespace lsm
  {
    template <
      typename TmplIterator
    >
    class IteratorPy
    {
    public:
      typedef TmplIterator                  Iterator;
      typedef typename Iterator::value_type value_type;

      IteratorPy(const Iterator &it);

      /**
       * Returns the next item in the sequence.
       */
      value_type next();

      static void exportIterator(
        const std::string &pythonName,
        const std::string &pythonDocReturnType="object"
      );

    private:
      Iterator m_it;
    };
  }
}

#include "Python/esys/lsm/geometry/IteratorPy.hpp"

#endif
