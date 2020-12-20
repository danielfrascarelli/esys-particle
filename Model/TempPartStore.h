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

#ifndef __TEMPPARTSTORE_H
#define __TEMPPARTSTORE_H

// --- Project includes ---
#include "vec3.h"

// --- STL includes ---
#include <map>
#include <utility>

using std::map;
using std::multimap;
using std::pair;
using std::make_pair;

/*!
  \class ATempPartStore
  \brief pure virtual base for TTempPartStore
*/
class ATempPartStore
{
 public:
  virtual ~ATempPartStore() {}
  virtual void addSlaveID(int,int,int,int)=0;
  virtual void addConnection(int,int,int)=0;
};

/*!
  \class TTempPartStore
  \brief class for the temporary storage and distribution of
  particle data

*/
template<typename T>
class TTempPartStore : public ATempPartStore
{
 private:
  multimap<int,T> m_mmap;
  map<int,typename multimap<int,T>::iterator> m_by_id;
  map<int,int> m_slave_id_map;

  double m_xmin,m_xsize,m_ymin,m_ysize,m_zmin,m_zsize;
  int m_nx,m_ny,m_nz;

  int coordsToIndex(int,int,int);
  int posToIndex(const Vec3&);

 public:
  TTempPartStore(const Vec3&,const Vec3&, int,int,int);

  virtual void addSlaveID(int,int,int,int);
  virtual void addParticle(const T&);
  virtual void addConnection(int,int,int);

  const multimap<int,T>& getMap() const {return m_mmap;};
};

#include "TempPartStore.hpp"

#endif //__TEMPPARTSTORE_H
