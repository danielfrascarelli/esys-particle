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

#ifndef __ROTELASTICINTERACTIONGROUP_H
#define __ROTELASTICINTERACTIONGROUP_H

#include "Model/InteractionGroup.h"
#include "Model/RotElasticInteraction.h"
#include "Model/IGParam.h"

//--- IO includes ---
#include <iostream>
using std::ostream;
using std::endl;


/*!
  \brief Class for a group of unbonded,elastic interactions
*/
template <class T>
class CRotElasticInteractionGroup : public APairInteractionGroup<T>
{
 protected:
  vector<CRotElasticInteraction> m_interactions;
  AParallelInteractionStorage* m_exIG; //<! if an interaction is in m_exIG, it can't be in m_interactions
  double m_kr; //<! Normal spring constant

 public:
  CRotElasticInteractionGroup();
  CRotElasticInteractionGroup(const CRotElasticIGP*);
  virtual ~CRotElasticInteractionGroup(){};

  virtual void setExIG(AParallelInteractionStorage* eg){m_exIG=eg;};
  void setParam(const CRotElasticIGP*);

  /**
   * Null op, don't require time step size.
   */
  virtual void setTimeStepSize(double dt)
  {
  }

  virtual void calcForces();
  virtual void Update(ParallelParticleArray<T>*);
  friend ostream& operator<< <>(ostream&,const CRotElasticInteractionGroup<T>&); 
};

#include "Model/RotElasticInteractionGroup.hpp"

#endif //__ELASTICINTERACTIONGROUP_H
