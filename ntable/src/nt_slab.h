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

#ifndef __NT_SLAB_H
#define __NT_SLAB_H

//--- project includes ---
#include "nts_iter.h"
#include "ntable.h"
#include "dslice.h"

//--- forward decls ---
template <typename T> 
class NTSlab;

template <typename T> bool operator==(const NTSlab<T>&,const NTSlab<T>&);
template <typename T> bool operator!=(const NTSlab<T>&,const NTSlab<T>&);
template <typename T> ostream& operator<<(ostream &, const NTSlab<T> &);

/*!
  \class NTSlab
  \brief representation of a slab of the search array of a NeigborTable
*/ 
template <typename T>
class NTSlab
{
 public: // types
  typedef NTSlab_iter<T> iterator;
  typedef NTSlab_iter<T> const_iterator; // dodgy
  typedef T value_type;

 private:
  NeighborTable<T> *m_table;
  DSlice m_sl;


 public:
  NTSlab();
  NTSlab(NeighborTable<T>*,DSlice);
  
  unsigned int slab_size()const {return m_sl.size();}
  unsigned int size()const;

  // begin and end iterators
  iterator begin();
  iterator end();

  iterator rbegin();
  iterator rend();

  //!< number of particles at a given gridpoint 
  unsigned int nparts_at_gridpoint(unsigned int idx) const {return m_table->nparts_at_gridpoint(m_sl[idx]);}; 
  
  // access ops
  T* ptr(typename NeighborTable<T>::indextype);
  T& ref(typename NeighborTable<T>::indextype);

  // insert ops
  void insert(iterator,const T&);
 
  // erase ops
  void erase(iterator);
  void erase(iterator,iterator);

  //comparison 
  friend bool operator== <>(const NTSlab&,const NTSlab&);
  friend bool operator!= <>(const NTSlab&,const NTSlab&);

  friend ostream& operator<< <>(ostream &, const NTSlab &);
};

#include "nt_slab.hpp"

#endif //__NT_SLAB_H
