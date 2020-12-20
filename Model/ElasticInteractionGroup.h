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

#ifndef __ELASTICINTERACTIONGROUP_H
#define __ELASTICINTERACTIONGROUP_H

#include "Model/InteractionGroup.h"
#include "Model/ElasticInteraction.h"
#include "Model/IGParam.h"

//--- IO includes ---
#include <iostream>
using std::ostream;
using std::endl;


/*!
  \brief Class for a group of unbonded,elastic interactions
*/
template <class T>
class CElasticInteractionGroup : public APairInteractionGroup<T>
{
 protected:
  vector<CElasticInteraction> m_interactions;
  AParallelInteractionStorage* m_exIG; //<! if an interaction is in m_exIG, it can't be in m_interactions
  double m_k; //<! spring constant
 
 public:
  CElasticInteractionGroup();
  CElasticInteractionGroup(const CElasticIGP*);
  virtual ~CElasticInteractionGroup(){};
  
  virtual void setExIG(AParallelInteractionStorage* eg){m_exIG=eg;};
  void setParam(const CElasticIGP*);

  virtual void calcForces();

  /**
   * Null operation, don't require time step size.
   */
  virtual void setTimeStepSize(double dt)
  {
  }
  
  virtual void Update(ParallelParticleArray<T>*);
  friend ostream& operator<< <>(ostream&,const CElasticInteractionGroup<T>&); 
};

#include "Model/ElasticInteractionGroup.hpp"

#endif //__ELASTICINTERACTIONGROUP_H
