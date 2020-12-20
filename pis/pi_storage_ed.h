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

#ifndef __PARALLEL_INTERACTION_STORAGE_ED_H
#define __PARALLEL_INTERACTION_STORAGE_ED_H

//--- project includes ---
#include "pi_storage_e.h"

//--- STL includes ---
#include <vector>
using std::vector;

//--- IO includes ---

/*!
  \class ParallelInteractionStorage_ED
  \brief parallel storage array with exchange for dynamically created interactions (friction etc.)
*/
template<typename P,typename I>
class ParallelInteractionStorage_ED : public ParallelInteractionStorage_E<P,I>
{
 public:
  //  typedef I ParallelInteractionStorage_ED::interaction_type;

 protected:
  int m_update_timestamp;
  vector<AParallelInteractionStorage*> m_exIG; //<! if an interaction is in m_exIG, it can't be in m_interactions

  bool isExcluded(const vector<int>); //<! check if a particle pair is in one of the excluding IGs 


 public:
  ParallelInteractionStorage_ED(AParallelParticleArray*,const typename I::ParameterType&);

  virtual void addExIG(AParallelInteractionStorage*);
  virtual bool update();
  virtual void setTimeStepSize(double dt);

  virtual void saveCheckPointData(std::ostream &oStream);
  virtual void loadCheckPointData(std::istream &iStream);

  virtual void calcHeatTrans();
  virtual void calcHeatFrict();

  virtual bool willSave(){ return true;};
};

#include "pis/pi_storage_ed.hpp"

#endif // __PARALLEL_INTERACTION_STORAGE_ED_H
