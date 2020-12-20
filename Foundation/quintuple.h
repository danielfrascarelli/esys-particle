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


#ifndef ESYS_LSMQUINTUPLE_H
#define ESYS_LSMQUINTUPLE_H

#include <boost/tuple/tuple.hpp>

namespace esys
{
  namespace lsm
  {
    template <typename T1, typename T2, typename T3, typename T4, typename T5>
    class quintuple : public boost::tuple<T1,T2,T3,T4,T5>
    {
    public:
      typedef boost::tuple<T1,T2,T3,T4,T5> inherited;
      inline quintuple() : inherited()
      {
      }

      inline quintuple(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5)
        : inherited(t1,t2,t3,t4,t5)
      {
      }

      inline quintuple(const quintuple &quin)
        : inherited(quin)
      {
      }

      inline quintuple &operator=(const quintuple &quin)
      {
        inherited::operator=(quin);
        return *this;
      }
    };
  }
}

#endif
