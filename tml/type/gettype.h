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

#ifndef __GETTYPE_H
#define __GETTYPE_H

//--- MPI ---
#include <mpi.h>

//--- STL includes ---
#include <utility>

#include <Foundation/triplet.h>
#include <Foundation/quadtuple.h>
#include <Foundation/quintuple.h>

using std::pair;

/*!
  \class SGetType
  \brief static function object to extract MPI type from data in a consistent way

*/
class SGetType
{
 public:
  template <typename T> MPI_Datatype operator()(const T &);

  template<typename T1,typename T2>
  MPI_Datatype operator()(const pair<T1,T2>&);

  template<typename T1, typename T2, typename T3>
  MPI_Datatype operator()(const esys::lsm::triplet<T1,T2,T3> &trip);

  template<typename T1, typename T2, typename T3, typename T4>
  MPI_Datatype operator()(const esys::lsm::quadtuple<T1,T2,T3,T4> &quad);

  template<typename T1, typename T2, typename T3, typename T4, typename T5>
  MPI_Datatype operator()(const esys::lsm::quintuple<T1,T2,T3,T4,T5> &quin);
};

/**
 * std::pair MPI-type initialisation stuff.
 */
template<typename T1,typename T2> 
struct tml_pair {
  static MPI_Datatype type;
  static bool initialized;
};
template<typename T1,typename T2> MPI_Datatype tml_pair<T1,T2>::type=MPI_DATATYPE_NULL;
template<typename T1,typename T2> bool tml_pair<T1,T2>::initialized=false;

/**
 * triplet MPI-type initialisation stuff.
 */
template<typename T1, typename T2, typename T3>
struct tml_trip {
  static MPI_Datatype type;
  static bool initialized;
};
template<typename T1, typename T2, typename T3> MPI_Datatype tml_trip<T1,T2,T3>::type=MPI_DATATYPE_NULL;
template<typename T1, typename T2, typename T3> bool tml_trip<T1,T2,T3>::initialized=false;

/**
 * Quadtuple MPI-type initialisation stuff.
 */
template<typename T1, typename T2, typename T3, typename T4>
struct tml_quad {
  static MPI_Datatype type;
  static bool initialized;
};
template<typename T1, typename T2, typename T3, typename T4> MPI_Datatype tml_quad<T1,T2,T3,T4>::type=MPI_DATATYPE_NULL;
template<typename T1, typename T2, typename T3, typename T4> bool tml_quad<T1,T2,T3,T4>::initialized=false;

/**
 * Quintuple MPI-type initialisation stuff.
 */
template<typename T1, typename T2, typename T3, typename T4, typename T5>
struct tml_quin {
  static MPI_Datatype type;
  static bool initialized;
};
template<typename T1, typename T2, typename T3, typename T4, typename T5> MPI_Datatype tml_quin<T1,T2,T3,T4,T5>::type=MPI_DATATYPE_NULL;
template<typename T1, typename T2, typename T3, typename T4, typename T5> bool tml_quin<T1,T2,T3,T4,T5>::initialized=false;

static SGetType GetType;

#include "tml/type/gettype.hpp"

#endif //__GETTYPE_H
