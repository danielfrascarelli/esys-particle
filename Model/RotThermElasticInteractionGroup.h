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

#ifndef __ROTTHERMELASTICINTERACTIONGROUP_H
#define __ROTTHERMELASTICINTERACTIONGROUP_H

#include "Foundation/console.h"
#include "InteractionGroup.h"
#include "RotThermElasticInteraction.h"
#include "IGParam.h"

//--- IO includes ---
#include <iostream>
using std::ostream;
using std::endl;


/*!
  Class for a group of unbonded,elastic interactions
*/
template <class T>
class CRotThermElasticInteractionGroup : public APairInteractionGroup<T>
{
 protected:
  vector<CRotThermElasticInteraction> m_interactions;
  AParallelInteractionStorage* m_exIG; //<! if an interaction is in m_exIG, it can't be in m_interactions
  double m_k; //<! spring constant
  double m_diffusivity  ;
 
 public:
  CRotThermElasticInteractionGroup();
  CRotThermElasticInteractionGroup(const CRotThermElasticIGP *);
  virtual ~CRotThermElasticInteractionGroup(){};
  
  virtual void setExIG(AParallelInteractionStorage* eg){m_exIG=eg;};
  void setParam(const CRotThermElasticIGP*);

  virtual void calcForces();
  virtual void Update(ParallelParticleArray<T>*);
  friend ostream& operator<< <>(ostream&,const CRotThermElasticInteractionGroup<T>&); 
};

#include "RotThermElasticInteractionGroup.hpp"

#endif //__ELASTICINTERACTIONGROUP_H
