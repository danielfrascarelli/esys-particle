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


#include <boost/python/iterator.hpp>

namespace esys
{
  namespace lsm
  {
    template<typename TmplIterator>
    IteratorPy<TmplIterator>::IteratorPy(const Iterator &it) : m_it(it)
    {
    }

    template<typename TmplIterator>
    typename IteratorPy<TmplIterator>::value_type
    IteratorPy<TmplIterator>::next()
    {
      if (m_it.hasNext())
      {
        return m_it.next();
      }
      boost::python::objects::stop_iteration_error();
      return m_it.next(); // to avoid compiler warning
    }

    template <typename TmplIterator>
    void IteratorPy<TmplIterator>::exportIterator(
      const std::string &pythonName,
      const std::string &pythonDocReturnType
    )
    {
      boost::python::class_<IteratorPy<TmplIterator> >(
        pythonName.c_str(),
        boost::python::no_init
      )
      .def(
#if PY_VERSION_HEX >= 0x03000000
        "__next__",
#else
        "next",
#endif
        &IteratorPy<TmplIterator>::next,
        boost::python::return_internal_reference<>(),
        (
          std::string("Returns the next object in the sequence.\n")
          +
          "@rtype: L{" + pythonDocReturnType + "}\n"
          +
          "@raise StopIteration: if there is no next element."
        ).c_str()
      )
      .def(
        "__iter__",
        boost::python::objects::identity_function(),
        "Returns self."
      )
      ;
    }
  }
}
