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

#ifndef __PARALLEL_INTERACTION_STORAGE_NE_T_H
#define __PARALLEL_INTERACTION_STORAGE_NE_T_H

//--- project includes ---
#include "pi_storage_ed.h"

//--- STL includes ---

//--- IO includes ---

/*!
  \class ParallelInteractionStorage_NE_T
  Parallel storage array with exchange for dynamically created interactions without exchange (elastic etc.) between particles with defined tags. An interaction is only created if particle tag & mask agree 
  

*/
template<typename P,typename I>
class ParallelInteractionStorage_NE_T : public ParallelInteractionStorage_NE<P,I>
{
 protected:
  int m_tag1, m_tag2; // particle tags
  int m_mask1, m_mask2; // tag masks
  
 public:
  ParallelInteractionStorage_NE_T(AParallelParticleArray*,const typename I::ParameterType&, int,int,int,int);

  virtual bool update();
  virtual bool willSave(){ return true;};
};

#include "pis/pi_storage_ne_t.hpp"

#endif // __PARALLEL_INTERACTION_STORAGE_NE_T_H
