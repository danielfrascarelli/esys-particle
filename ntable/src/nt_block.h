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

#ifndef __NT_BLOCK_H
#define __NT_BLOCK_H

//--- project includes ---
#include "ntable/src/ntb_iter.h"
#include "ntable/src/ntable.h"

//--- forward decls ---
template <typename T> 
class NTBlock;

template <typename T> bool operator==(const NTBlock<T>&,const NTBlock<T>&);
template <typename T> bool operator!=(const NTBlock<T>&,const NTBlock<T>&);
template <typename T> ostream& operator<<(ostream &, const NTBlock<T>&);
/*!
  \class NTBlock
  \brief representation of a slab of the search array of a NeigborTable
*/ 
template <typename T>
class NTBlock
{
 public: // types
  typedef NTBlock_iter<T> iterator;
  typedef NTBlock_iter<T> const_iterator; // dodgy
  typedef T value_type;

 private:
  NeighborTable<T> *m_table;
  int m_xmin,m_xmax,m_ymin,m_ymax,m_zmin,m_zmax;

 public:
  NTBlock();
  NTBlock(NeighborTable<T>*,int,int,int,int,int,int);
  
  unsigned int size();

  // begin and end iterators
  iterator begin();
  iterator end();

  /* iterator rbegin(); */
/*   iterator rend(); */

  //!< number of particles at a given gridpoint 
  unsigned int nparts_at_gridpoint(int x,int y,int z) const {return m_table->nparts_at_gridpoint(m_table->index(x,y,z));}; 
  
  // access ops
  T* ptr(int,int,int,int);
  T& ref(int,int,int,int);

  //comparison 
  friend bool operator== <>(const NTBlock&,const NTBlock&);
  friend bool operator!= <>(const NTBlock&,const NTBlock&);

  friend ostream& operator<< <>(ostream &, const NTBlock &);

  friend class NTBlock_iter<T>;
};

#include "ntable/src/nt_block.hpp"

#endif //__NT_BLOCK_H
