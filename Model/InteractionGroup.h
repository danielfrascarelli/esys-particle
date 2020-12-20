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

#ifndef __INTERACTIONGROUP_H
#define __INTERACTIONGROUP_H

// -- project includes -- 
#include "Model/IGParam.h"
#include "Model/InteractionParam.h"

// -- STL includes --
#include <set>
#include <utility>

using std::set;
using std::pair;

class AParallelInteractionStorage;
template <class T> class ParallelParticleArray;

/*!
  \brief Abstract base class for a group of interactions
*/
template <class T>
class AInteractionGroup
{
public:
  virtual ~AInteractionGroup(){};

  virtual void Update(ParallelParticleArray<T>*)=0;
  virtual void calcForces()=0;

  virtual void setTimeStepSize(double dt) = 0;
};

/*!
  \brief Abstract base class for a group of pair interactions

  The difference to AInteractionGroup is the existence of a function 
  bool isIn(int,int) which returns if an interaction between particles
  with the given Ids is in this group. (Necessary because bonded -> not 
  elastic,frictional etc.)
*/
template<class T>
class APairInteractionGroup : public AInteractionGroup<T>
{
 protected:
  /*!
    having an extra set of all pairs if particle-Ids of the interactions
    in the group costs some memory, but speeds up isIn to O(log N). It would 
    be O(N) if implemented with a find() over the vector of inteactions.
    \todo replace set (O(logN)) with hashset (O(1))
   */
  set<pair<int,int> > m_set;
  unsigned int m_update_timestamp;

 public: 
  bool isIn(int,int);
  virtual void setExIG(AParallelInteractionStorage* eg){};
};

#include "Model/InteractionGroup.hpp"

#endif //__INTERACTIONGROUP_H
