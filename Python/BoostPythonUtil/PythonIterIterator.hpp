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


namespace esys
{
  namespace lsm
  {
    namespace bpu
    {
      template <typename TmplExtractType>
      PythonIterIterator<TmplExtractType>::PythonIterIterator(
        boost::python::object &iteratable
      )
        : m_hasNext(true),
          m_next(),
          m_iter(iteratable.attr("__iter__")())
      {
        update();
      }

      template <typename TmplExtractType>
      bool PythonIterIterator<TmplExtractType>::hasNext() const
      {
        return m_hasNext;
      }

      template <typename TmplExtractType>
      TmplExtractType PythonIterIterator<TmplExtractType>::next()
      {
        boost::python::object next = m_next;
        update();
        return boost::python::extract<TmplExtractType>(next);
      }

      template <typename TmplExtractType>
      void PythonIterIterator<TmplExtractType>::update()
      {
        try
        {
#if PY_VERSION_HEX >= 0x03000000
          m_next = m_iter.attr("__next__")();
#else
          m_next = m_iter.attr("next")();
#endif
        }
        catch (boost::python::error_already_set &e)
        {
          if (!PyErr_ExceptionMatches(PyExc_StopIteration))
          {
            throw;
          }
          PyErr_Clear();
          m_hasNext = false;
        }
      }
    }
  }
}
