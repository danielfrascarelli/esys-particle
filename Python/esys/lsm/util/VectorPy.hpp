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
    template <typename TE>
    VectorPy<TE>::VectorPy() : Inherited()
    {
    }

    template <typename TE>
    VectorPy<TE>::VectorPy(const VectorPy &vec) : Inherited(vec)
    {
    }

    template <typename TE>
    VectorPy<TE>::VectorPy(const Inherited &vec) : Inherited(vec)
    {
    }

    template <typename TE>
    VectorPy<TE>::VectorPy(boost::python::object &iterable)
      : Inherited()
    {
      bpu::PythonIterIterator<TE> it(iterable);
      while (it.hasNext())
      {
        this->push_back(it.next());
      }
    }

    template <typename TE>
    void VectorPy<TE>::append(const_reference e)
    {
      this->push_back(e);
    }

    template <typename TE>
    size_t VectorPy<TE>::getIndex(int i) const
    {
      if (i < 0)
      {
          i += this->size();
      }
      if ((i >= static_cast<int>(this->size())) || (i < 0))
      {
        PyErr_SetString(PyExc_IndexError, "Index out of range");
        boost::python::throw_error_already_set();
      }
      return static_cast<size_t>(i);
    }

    template <typename TE>
    typename VectorPy<TE>::reference VectorPy<TE>::getItem(int i)
    {
      return (*this)[getIndex(i)];
    }

    template <typename TE>
    void VectorPy<TE>::setItem(int i, const_reference e)
    {
      (*this)[getIndex(i)] = e;
    }

    template <typename TE>
    boost::python::tuple
    VectorPy<TE>::PickleSuite::getinitargs(VectorPy const& v)
    {
      return boost::python::make_tuple(boost::python::list(v));
    }

    template <typename TE>
    boost::python::class_<VectorPy<TE> >
    VectorPy<TE>::exportVector(
      const std::string &pyClassName,
      const std::string &pyClassDocString
    )
    {
      return
        boost::python::class_<VectorPy<TE> >(
          pyClassName.c_str(),
          pyClassDocString.c_str(),
          boost::python::init<>()
        )
        .def(
          boost::python::init<const VectorPy<TE> &>(
            (boost::python::arg("iterable"))
          )
        )
        .def(
          boost::python::init<boost::python::object &>(
            (boost::python::arg("iterable")),
            (
              std::string() + 
              "Constructs vec of elements from specifed sequence.\n" +
              "@type iterable: iterable\n" +
              "@kwarg iterable: copy and insert elements from C{iterable}" +
              " into this vec.\n"
            ).c_str()
          )
        )
        .def(
          "__iter__",
          boost::python::iterator<
            VectorPy,
            boost::python::return_internal_reference<>
          >()
        )
        .def(
          "__len__",
          &VectorPy<TE>::size
        )
        .def(
          "__getitem__",
          &VectorPy<TE>::getItem,
          boost::python::return_internal_reference<>()
        )
        .def(
          "__setitem__",
          &VectorPy<TE>::setItem
        )
        .def(
          "append",
          &VectorPy<TE>::append,
          (boost::python::arg("elem")),
          "Appends an element to the end of this vector.\n"
          "@type elem: object\n"
          "@kwarg elem: Copy this element as the new last element."
        )
        .def(
          "clear",
          &VectorPy<TE>::clear,
          "Removes all elements from this vector."
        )
        .def(
          boost::python::self == boost::python::self
        )
        .def_pickle(typename VectorPy::PickleSuite())
        ;
    }
  }
}
