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

#ifndef __PARALLEL_INTERACTION_STORAGE_NE_H
#define __PARALLEL_INTERACTION_STORAGE_NE_H

//--- project includes ---
#include "pis/pi_storage.h"

//--- STL includes ---
#include <vector>

//--- IO includes ---

/*!
  \brief parallel storage array without exchange for dynamically created interactions (elastic)
*/
template<typename P,typename I>
class ParallelInteractionStorage_NE : public TParallelInteractionStorage<I>
{
 protected:
  int m_update_timestamp;
  vector<AParallelInteractionStorage*> m_exIG; //<! if an interaction is in m_exIG, it can't be in m_interactions

  set<pair<int,int> > m_set; // evil hack, should be vector<int>, not pair<int,int>
  typename I::ParameterType m_param;

  bool isExcluded(const vector<int>); //<! check if a particle pair is in one of the excluding IGs 

 public:
  ParallelInteractionStorage_NE(AParallelParticleArray*,const typename I::ParameterType&);

  virtual void addExIG(AParallelInteractionStorage*);
  virtual bool update();
  virtual void exchange(){}; //!< do nothing
  virtual void rebuild(){}; //!< do nothing
  virtual void tryInsert(const vector<int>&){};//!< do nothing
  virtual bool isIn(const vector<int>&);
  virtual void calcForces();
  virtual void setTimeStepSize(double dt)
  {
  }
  virtual void calcHeatTrans();
};

#include "pis/pi_storage_ne.hpp"

#endif // __PARALLEL_INTERACTION_STORAGE_NE_H
