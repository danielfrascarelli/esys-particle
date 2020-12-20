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

#ifndef ESYS_LSMVECTORPY_H
#define ESYS_LSMVECTORPY_H

#include <boost/python.hpp>
#include <vector>
#include <iterator>

namespace esys
{
  namespace lsm
  {
    template <typename TmplElem>
    class VectorPy : public std::vector<TmplElem>
    {
    public:
      typedef std::vector<TmplElem> Inherited;
      typedef typename Inherited::const_reference const_reference;
      typedef typename Inherited::reference reference;

      class PickleSuite : public boost::python::pickle_suite
      {
      public:
        static
        boost::python::tuple
        getinitargs(VectorPy const& v);
      };

      VectorPy();

      VectorPy(const VectorPy &vec);

      VectorPy(const Inherited &vec);

      VectorPy(boost::python::object &iterable);

      void append(const_reference e);

      size_t getIndex(int i) const;

      reference getItem(int i);

      void setItem(int i, const_reference e);

      static boost::python::class_<VectorPy>
        exportVector(
          const std::string &pyClassName,
          const std::string &pyClassDocString
        );
    };
  }
}

#include "Python/esys/lsm/util/VectorPy.hpp"

#endif
