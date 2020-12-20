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

#ifndef ESYS_LSMSETPY_H
#define ESYS_LSMSETPY_H

#include <boost/python.hpp>
#include <set>
#include <iterator>

namespace esys
{
  namespace lsm
  {
    template <typename TmplElem, typename TmplCompare=std::less<TmplElem> >
    class SetPy : public std::set<TmplElem, TmplCompare>
    {
    public:
      typedef std::set<TmplElem, TmplCompare> Inherited;

      class PickleSuite : public boost::python::pickle_suite
      {
      public:
        static
        boost::python::tuple
        getinitargs(SetPy const& s);
      };

      SetPy();

      SetPy(const SetPy &set);

      SetPy(const Inherited &set);

      SetPy(boost::python::object &iterable);

      SetPy getUnion(const SetPy &set) const;

      SetPy getDifference(const SetPy &set) const;

      SetPy getIntersection(const SetPy &set) const;

      static boost::python::class_<SetPy>
        exportSet(
          const std::string &pyClassName,
          const std::string &pyClassDocString
        );

      protected:
        typedef std::insert_iterator<SetPy> InsertIterator;
    };
  }
}

#include "Python/esys/lsm/util/SetPy.hpp"

#endif
