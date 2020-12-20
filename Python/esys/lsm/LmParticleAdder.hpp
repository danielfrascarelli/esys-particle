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


#include <boost/mpl/placeholders.hpp>
#include "Python/BoostPythonUtil/Util.h"
#include "Python/BoostPythonUtil/PythonIterIterator.h"
#include <stdexcept>

using namespace boost::mpl::placeholders;
namespace esys
{
  namespace lsm
  {
    template <class T>
    struct Wrap
    {
    };
    
    template <class WrappedT>
    class PtrWrap
    {
    public:
      PtrWrap(WrappedT &t) : m_p(&t)
      {
      }

      template <typename T>
      void operator()(Wrap<T> wrappedT)
      {
        (*m_p)(wrappedT);
      }
    private:
      WrappedT *m_p;
    };

    class ExtractIndexer
    {
    public:
      ExtractIndexer(boost::python::object pyOb)
        : m_currIndex(-1),
          m_extractIndex(-1),
          m_pyOb(pyOb)
      {
      }

      int getExtractIndex() const
      {
        return m_extractIndex;
      }

      template <class T>
      void operator()(Wrap<T>)
      {
        ++m_currIndex;
        if (m_extractIndex < 0)
        {
          boost::python::extract<T &> extractor(m_pyOb);
          if (extractor.check())
          {
            m_extractIndex = m_currIndex;
          }
        }
      }
      
    private:
      int m_currIndex;
      int m_extractIndex;
      boost::python::object m_pyOb;
    };
    
    typedef bpu::PythonIterIterator<boost::python::object> PyObjectIterator;
    typedef std::vector<boost::python::list> PyListVector;

    template <class TmplLsmParticle>
    class LmAdder
    {
    public:
      LmAdder(CLatticeMaster &lm, PyListVector &pyListVector)
        : m_currIndex(-1),
          m_pLm(&lm),
          m_pPyListVector(&pyListVector)
      {
      }

      template <class T>
      void operator()(Wrap<T>)
      {
        ++(this->m_currIndex);

        typedef bpu::PythonIterIterator<T &> TRefIterator;
        TRefIterator it((*m_pPyListVector)[m_currIndex]);
        this->m_pLm->template addParticles<TRefIterator,TmplLsmParticle>(it);
      }

    private:
      int            m_currIndex;
      CLatticeMaster *m_pLm;
      PyListVector   *m_pPyListVector;
    };
    
    template <class TmplMplVector, class TmplLsmParticle>
    void LmParticleAdder<TmplMplVector,TmplLsmParticle>::addParticles(
      boost::python::object &iterable,
      CLatticeMaster &lm
    )
    {
      PyObjectIterator it(iterable);
      const size_t numTypes = boost::mpl::size<MplVector>::type::value;
      PyListVector pyListVector;
      pyListVector.reserve(numTypes);
      for (size_t i = 0; i < numTypes; i++)
      {
        pyListVector.push_back(boost::python::list());
      }
      while (it.hasNext())
      {
        boost::python::object pyOb = it.next();
        ExtractIndexer extractIndex(pyOb);
        boost::mpl::for_each<MplVector, Wrap<boost::mpl::placeholders::_1> >(
          PtrWrap<ExtractIndexer>(extractIndex)
        );
        if (extractIndex.getExtractIndex() >= 0)
        {
          pyListVector[extractIndex.getExtractIndex()].append(pyOb);
        }
        else
        {
          throw std::runtime_error(
            boost::python::extract<std::string>(
              std::string("Could not extract C++ type from python object:")
              +
              boost::python::str(pyOb)
            )
          );
        }
      }
      boost::mpl::for_each<MplVector, Wrap<boost::mpl::placeholders::_1> >(LmAdder<LsmParticle>(lm, pyListVector));
    }
  }
}
