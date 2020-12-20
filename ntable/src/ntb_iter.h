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

#ifndef __NTB_BLOCK_H
#define __NTB_BLOCK_H

//--- project includes ---
#include "ntable/src/ntable.h"

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
class NTBlock_iter;

template <typename T> bool operator==(const NTBlock_iter<T>&,const NTBlock_iter<T>&);
template <typename T> bool operator!=(const NTBlock_iter<T>&,const NTBlock_iter<T>&);

template <typename T> 
class NTBlock;

/*!
  \class NTBlock_iter
  \brief iterator for a NTBlock
*/
template <typename T>
class NTBlock_iter
{
 private:
  NTBlock<T> *m_block;
  int m_ix,m_iy,m_iz,m_ig;

 public:
  NTBlock_iter(NTBlock<T>*,int,int,int,int);

  //! move ops
  NTBlock_iter& operator++();
  NTBlock_iter operator++(int);

  //! access ops
  T* operator->();
  T& operator*();

  //comparison 
  friend bool operator== <>(const NTBlock_iter&,const NTBlock_iter&);
  friend bool operator!= <>(const NTBlock_iter&,const NTBlock_iter&);
};

#include "ntable/src/ntb_iter.hpp"

#endif //__NTB_ITER_H

