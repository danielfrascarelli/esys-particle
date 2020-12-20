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

#ifndef __PARALLEL_INTERACTION_STORAGE_E_H
#define __PARALLEL_INTERACTION_STORAGE_E_H

//--- project includes ---
#include "pi_storage.h"
#include "tml/comm/cart_comm.h"

//--- STL includes ---
#include <utility>
#include <set>
#include <vector>
using std::vector;
using std::pair;
using std::make_pair;
using std::set;

class AParallelParticleArray;

/*!
  \class ParallelInteractionStorage_E
  \brief parallel interaction storage array with exchange
*/
template<typename P,typename I>
class ParallelInteractionStorage_E : public TParallelInteractionStorage<I>
{
 public: // types
  //  typedef I ParallelInteractionStorage_E::interaction_type;
  typedef TParallelInteractionStorage<I>          Inherited;
  typedef typename Inherited::InteractionIterator InteractionIterator;
  bool m_unbreakable;

 private:

  static const int m_exchg_tag;
  
  void exchange_boundary(int,int);

 protected:
  TML_CartComm m_comm;
  set<pair<int,int> > m_set; // evil hack, should be std::vector<int>, not pair<int,int>
  typename I::ParameterType m_param;

 public:
  ParallelInteractionStorage_E(AParallelParticleArray *, const typename I::ParameterType &);

  virtual void setUnbreakable(bool);
  virtual void exchange();
  virtual void rebuild();
  virtual void tryInsert(const I&);
  virtual void tryInsert(const std::vector<int>&);
  virtual bool isIn(const std::vector<int>&);
  virtual void setTimeStepSize(double){} //!< does nothing
  virtual void calcForces();
};

#include "pis/pi_storage_e.hpp"

#endif // __PARALLEL_INTERACTION_STORAGE_E_H
