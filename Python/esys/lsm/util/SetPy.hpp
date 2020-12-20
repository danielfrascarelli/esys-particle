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

#include "Python/BoostPythonUtil/PythonIterIterator.h"
#include <algorithm>

namespace esys
{
  namespace lsm
  {
    template <typename TE, typename TC>
    SetPy<TE,TC>::SetPy() : Inherited()
    {
    }

    template <typename TE, typename TC>
    SetPy<TE,TC>::SetPy(const SetPy &set) : Inherited(set)
    {
    }

    template <typename TE, typename TC>
    SetPy<TE,TC>::SetPy(const Inherited &set) : Inherited(set)
    {
    }

    template <typename TE, typename TC>
    SetPy<TE,TC>::SetPy(boost::python::object &iterable)
      : Inherited()
    {
      bpu::PythonIterIterator<TE> it(iterable);
      while (it.hasNext())
      {
        this->insert(it.next());
      }
    }

    template <typename TE, typename TC>
    SetPy<TE,TC> SetPy<TE,TC>::getUnion(const SetPy &set) const
    {
      SetPy resultSet;
      std::set_union(
        this->begin(),
        this->end(),
        set.begin(),
        set.end(),
        InsertIterator(resultSet, resultSet.begin())
      );
      return resultSet;
    }

    template <typename TE, typename TC>
    SetPy<TE,TC> SetPy<TE,TC>::getDifference(const SetPy &set) const
    {
      SetPy resultSet;
      std::set_difference(
        this->begin(),
        this->end(),
        set.begin(),
        set.end(),
        InsertIterator(resultSet, resultSet.begin())
      );
      return resultSet;   
    }

    template <typename TE, typename TC>
    SetPy<TE,TC> SetPy<TE,TC>::getIntersection(const SetPy &set) const
    {
      SetPy resultSet;
      std::set_intersection(
        this->begin(),
        this->end(),
        set.begin(),
        set.end(),
        InsertIterator(resultSet, resultSet.begin())
      );
      return resultSet;
    }

    template <typename TE, typename TC>
    boost::python::tuple
    SetPy<TE,TC>::PickleSuite::getinitargs(SetPy const& s)
    {
      return boost::python::make_tuple(boost::python::list(s));
    }

    template <typename TE, typename TC>
    boost::python::class_<SetPy<TE,TC> >
    SetPy<TE,TC>::exportSet(
      const std::string &pyClassName,
      const std::string &pyClassDocString
    )
    {
      return
        boost::python::class_<SetPy<TE,TC> >(
          pyClassName.c_str(),
          pyClassDocString.c_str(),
          boost::python::init<>()
        )
        .def(
          boost::python::init<const SetPy<TE,TC> &>(
            (boost::python::arg("iterable"))
          )
        )
        .def(
          boost::python::init<boost::python::object &>(
            (boost::python::arg("iterable")),
            (
              std::string() + 
              "Constructs set of elements from specifed sequence.\n" +
              "@type iterable: iterable\n" +
              "@kwarg iterable: copy and insert elements from C{iterable}" +
              " into this set.\n"
            ).c_str()
          )
        )
        .def(
          "union",
          &SetPy::getUnion,
          (boost::python::arg("set")),
          (
            std::string() +
            "Returns new set which is the union of this set and C{set}.\n" +
            "@type set: " + pyClassName + "\n" +
            "@kwarg set: Set with which union is performed.\n" +
            "@rtype: " + pyClassName + "\n"
            "@return: new set which is the union of two sets."
          ).c_str()
        )
        .def(
          "difference",
          &SetPy::getDifference,
          (boost::python::arg("set")),
          (
            std::string() +
            "Returns new set which contains elements of this set not " + 
            " contained in the specified C{set}.\n" +
            "@type set: " + pyClassName + "\n" +
            "@kwarg set: Set with which difference is performed.\n" +
            "@rtype: " + pyClassName + "\n"
            "@return: new set which is the difference of two sets."
          ).c_str()
        )
        .def(
          "intersection",
          &SetPy::getIntersection,
          (boost::python::arg("set")),
          (
            std::string() +
            "Returns new set which contains elements common to this set " +
            " and the specified C{set}.\n" +
            "@type set: " + pyClassName + "\n" +
            "@kwarg set: Set with which intersection is performed.\n" +
            "@rtype: " + pyClassName + "\n"
            "@return: new set which is the intersection of two sets."
          ).c_str()
        )
        .def(
          "__iter__",
          boost::python::iterator<
            SetPy,
            boost::python::return_internal_reference<>
          >()
        )
        .def(
          "__len__",
          &SetPy::size
        )
        .def(
          boost::python::self == boost::python::self
        )
        .def_pickle(typename SetPy::PickleSuite())
        ;
    }
  }
}
