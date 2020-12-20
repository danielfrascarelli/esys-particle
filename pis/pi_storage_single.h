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

#ifndef __PARALLEL_INTERACTION_STORAGE_SINGLE_H
#define __PARALLEL_INTERACTION_STORAGE_SINGLE_H

#include "pis/pi_storage.h"
#include <vector>

//--- IO includes ---

/*!
  \brief parallel storage array without exchange for dynamically created single particle
  interactions (i.e. Damping...)
*/
template<typename P,typename I>
class ParallelInteractionStorage_Single : public TParallelInteractionStorage<I>
{
 protected:
  typename I::ParameterType m_param;

 public:
  ParallelInteractionStorage_Single(AParallelParticleArray*,const typename I::ParameterType&);

  virtual void addExIG(AParallelInteractionStorage*){}; // do nothing
  virtual bool update();
  virtual void exchange(){}; //!< do nothing
  virtual void rebuild(){}; //!< do nothing
  virtual void tryInsert(const vector<int>&){};//!< do nothing
  virtual bool isIn(const vector<int>&){return true;}; 
  virtual void calcForces();
  virtual void setTimeStepSize(double dt);
};

#include "pi_storage_single.hpp"

#endif //__PARALLEL_INTERACTION_STORAGE_SINGLE_H
