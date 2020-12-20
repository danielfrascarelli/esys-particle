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

#ifndef __PARALLEL_INTERACTION_STORAGE_EB_H
#define __PARALLEL_INTERACTION_STORAGE_EB_H

//--- project includes ---
#include "pis/pi_storage_e.h"

//--- STL includes ---
#include <vector>

//--- IO includes ---

/*!
  \class ParallelInteractionStorage_EB
  \brief parallel storage array with exchange for bonded/breakable interactions 
*/
template<typename P,typename I>
class ParallelInteractionStorage_EB : public ParallelInteractionStorage_E<P,I>
{
 public:
  //  typedef I ParallelInteractionStorage_EB::interaction_type;
  typedef ParallelInteractionStorage_E<P,I>       Inherited;
  typedef typename Inherited::InteractionIterator InteractionIterator;
  bool m_unbreakable;

 public:
  ParallelInteractionStorage_EB(AParallelParticleArray*,const typename I::ParameterType&);

  virtual bool update();
  void setUnbreakable(bool);
  virtual void calcHeatTrans();
  virtual void saveCheckPointData(std::ostream &oStream);
  virtual void loadCheckPointData(std::istream &iStream);
  virtual void saveSnapShotData(std::ostream&);
};

#include "pis/pi_storage_eb.hpp"

#endif // __PARALLEL_INTERACTION_STORAGE_EB_H
