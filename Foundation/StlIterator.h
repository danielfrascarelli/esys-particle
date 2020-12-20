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


#ifndef ESYS_LSMSTLITERATOR_H
#define ESYS_LSMSTLITERATOR_H

#include <string>
#include <vector>
#include <sstream>

namespace esys
{
  namespace lsm
  {
    template <typename TmplIteratable>
    class ForwardConstIterator;

    template <typename TmplIteratable>
    class ForwardIterator
    {
    public:
      typedef typename TmplIteratable::value_type value_type;
      typedef typename TmplIteratable::reference  reference;
      typedef typename TmplIteratable::iterator   iterator;

      inline ForwardIterator(TmplIteratable &container)
        : m_it(container.begin()),
          m_end(container.end())
      {
      }

      inline ForwardIterator(const iterator &begin, const iterator &end)
        : m_it(begin),
          m_end(end)
      {
      }

      inline bool hasNext() const
      {
        return (m_it != m_end);
      }

      inline reference current() const
      {
        iterator it = m_it;
        it--;
        return *(it);
      }

      inline reference next()
      {
        reference ref = (*m_it);
        m_it++;
        return ref;
      }

      friend class ForwardConstIterator<TmplIteratable>;

    private:
      iterator m_it;
      iterator m_end;
    };
    
    template <typename TmplIteratable>
    class ForwardConstIterator
    {
    public:
      typedef typename TmplIteratable::value_type       value_type;
      typedef typename TmplIteratable::const_reference  reference;
      typedef typename TmplIteratable::const_iterator   iterator;

      inline ForwardConstIterator(const TmplIteratable &container)
        : m_it(container.begin()),
          m_end(container.end())
      {
      }

      inline ForwardConstIterator(const iterator &begin, const iterator &end)
        : m_it(begin),
          m_end(end)
      {
      }

      inline ForwardConstIterator(const ForwardIterator<TmplIteratable> &it)
        : m_it(it.m_it),
          m_end(it.m_end)
      {
      }

      inline bool hasNext() const
      {
        return (m_it != m_end);
      }

      inline reference current() const
      {
        iterator it = m_it;
        it--;
        return *(it);
      }

      inline reference next()
      {
        reference ref = (*m_it);
        m_it++;
        return ref;
      }

    private:
      iterator m_it;
      iterator m_end;
    };
  }
}

#endif
