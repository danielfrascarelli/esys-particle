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

#ifndef __NTS_ITER_H
#define __NTS_ITER_H

//--- project includes ---
#include "ntable.h"

//---system includes---
#include <cstddef>

//--- STL includes ---
#include <valarray>
#include <vector>
#include <list>

using std::slice;
using std::vector;
using std::list;

//--- forward decls ---
template <typename T> 
class NTSlab_iter;

template <typename T> bool operator==(const NTSlab_iter<T>&,const NTSlab_iter<T>&);
template <typename T> bool operator!=(const NTSlab_iter<T>&,const NTSlab_iter<T>&);

template <typename T> 
class NTSlab;

/*!
  \class NTSlab_iter
  \brief iterator for a NTSlab
*/
template <typename T>
class NTSlab_iter
{
 private:
  NTSlab<T> *m_slab;
  typename NeighborTable<T>::indextype m_curr;

 public:
  NTSlab_iter(NTSlab<T>*,typename NeighborTable<T>::indextype);

  //! move ops
  NTSlab_iter& operator++();
  NTSlab_iter operator++(int);

  NTSlab_iter& operator--();
  NTSlab_iter operator--(int);


  //! access ops
  T* operator->();
  T& operator*();

  typename NeighborTable<T>::indextype index() const;

  //comparison 
  friend bool operator== <>(const NTSlab_iter&,const NTSlab_iter&);
  friend bool operator!= <>(const NTSlab_iter&,const NTSlab_iter&);
};

#include "nts_iter.hpp"

#endif //__NTS_ITER_H
