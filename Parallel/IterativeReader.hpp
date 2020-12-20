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
    template <class TmplData>
    IStreamIterator<TmplData>::IStreamIterator(std::istream &iStream, int numParticles)
      : m_pIStream(&iStream),
        m_numRemaining(numParticles)
    {
    }

    template <class TmplData>
    IStreamIterator<TmplData>::~IStreamIterator()
    {
    }

    template <class TmplData>
    bool IStreamIterator<TmplData>::hasNext() const
    {
      return (m_numRemaining > 0);
    }

    template <class TmplData>
    void IStreamIterator<TmplData>::readDataFromStream()
    {
      m_data.read(*m_pIStream);
    }

    template <class TmplData>
    const TmplData &IStreamIterator<TmplData>::next()
    {
      m_numRemaining--;
      readDataFromStream();
      return m_data;
    }

    template <class TmplData>
    int IStreamIterator<TmplData>::getNumRemaining() const
    {
      return m_numRemaining;
    }

//===============================================================================
    template <class TmplIterator>
    IterativeReader<TmplIterator>::IterativeReader(std::istream &iStream)
      : m_numElements(-1),
        m_pIStream(&iStream),
        m_iteratorPtr()
    {
    }

    template <class TmplIterator>
    IterativeReader<TmplIterator>::~IterativeReader()
    {
    }
        
    template <class TmplIterator>
    int IterativeReader<TmplIterator>::getNumElements() const
    {
      return m_numElements;
    }

    template <class TmplIterator>
    void IterativeReader<TmplIterator>::setNumElements(int numElements)
    {
      m_numElements = numElements;
    }

    template <class TmplIterator>
    std::istream &IterativeReader<TmplIterator>::getIStream()
    {
      return *m_pIStream;
    }

    template <class TmplIterator>
    const std::istream &IterativeReader<TmplIterator>::getIStream() const
    {
      return *m_pIStream;
    }

    template <class TmplIterator>
    TmplIterator *IterativeReader<TmplIterator>::createNewIterator()
    {
      return new TmplIterator(*m_pIStream, getNumElements());
    }

    template <class TmplIterator>
    void IterativeReader<TmplIterator>::initialise()
    {
      m_iteratorPtr = IteratorAutoPtr(createNewIterator());
    }

    template <class TmplIterator>
    bool IterativeReader<TmplIterator>::isInitialised() const
    {
      return (m_iteratorPtr.get() != NULL);
    }
    
    template <class TmplIterator>
    typename IterativeReader<TmplIterator>::Iterator &IterativeReader<TmplIterator>::getIterator()
    {
      if (!isInitialised())
      {
        initialise();
      }
      return (*(m_iteratorPtr));
    }
  }
}
