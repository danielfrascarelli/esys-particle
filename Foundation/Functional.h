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

#ifndef ESYS_LSMFUNCTIONAL_H
#define ESYS_LSMFUNCTIONAL_H

#if HAVE_CONFIG_H
#include "config.h"
#endif

#if HAVE_EXT_FUNCTIONAL
#include <ext/functional>
#endif
#if HAVE_FUNCTIONAL
#include <functional>
#endif

#if HAVE_STD__SELECT1ST_PAIR_
namespace ext {
template <class _Pair> struct select1st
  : public std::select1st<_Pair> {};
}
#elif HAVE___GNU_CXX__SELECT1ST_PAIR_
namespace ext {
template <class _Pair> struct select1st
  : public __gnu_cxx::select1st<_Pair> {};
}
#elif !HAVE_EXT__SELECT1ST_PAIR_
namespace ext
{
  template <class _Pair>
  struct select1st
    : public std::unary_function<_Pair, typename _Pair::first_type>
  {
    typename _Pair::first_type&
    operator()(_Pair& __x) const
    { return __x.first; }

    const typename _Pair::first_type&
    operator()(const _Pair& __x) const
    { return __x.first; }
  };
}
#endif

#if HAVE_STD__SELECT2ND_PAIR_
namespace ext {
template <class _Pair> struct select2nd
  : public std::select2nd<_Pair> {};
}
#elif HAVE___GNU_CXX__SELECT2ND_PAIR_
namespace ext {
template <class _Pair> struct select2nd
  : public __gnu_cxx::select2nd<_Pair> {};
}
#elif !HAVE_EXT__SELECT2ND_PAIR_
namespace ext
{
  template <class _Pair>
  struct select2nd
    : public std::unary_function<_Pair, typename _Pair::second_type>
  {
    typename _Pair::second_type&
    operator()(_Pair& __x) const
    { return __x.second; }

    const typename _Pair::second_type&
    operator()(const _Pair& __x) const
    { return __x.second; }
  };
}
#endif

#endif
