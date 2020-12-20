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


#ifndef ESYS_LSMQUADTUPLE_H
#define ESYS_LSMQUADTUPLE_H

#include <boost/tuple/tuple.hpp>

namespace esys
{
  namespace lsm
  {
    template <typename T1, typename T2, typename T3, typename T4>
    class quadtuple : public boost::tuple<T1,T2,T3,T4>
    {
    public:
      typedef boost::tuple<T1,T2,T3,T4> inherited;
      inline quadtuple() : inherited()
      {
      }

      inline quadtuple(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4)
        : inherited(t1,t2,t3,t4)
      {
      }

      inline quadtuple(const quadtuple &quad)
        : inherited(quad)
      {
      }

      inline quadtuple &operator=(const quadtuple &quad)
      {
        inherited::operator=(quad);
        return *this;
      }
    };
  }
}

#endif
