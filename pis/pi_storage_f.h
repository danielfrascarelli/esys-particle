/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2014 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0          //
//                                                         //
/////////////////////////////////////////////////////////////


#ifndef __PARALLEL_INTERACTION_STORAGE_F_H
#define __PARALLEL_INTERACTION_STORAGE_F_H

#include "Foundation/vec3.h"
#include "ppa/src/pp_array.h"
#include "Model/FluidInteraction.h"
#include "pis/pi_storage.h"

#include <list>

using std::list;

class AParallelParticleArray;


/*!
  \class ParallelInteractionStorage_F
  \brief class for parallel fluid interaction storage array.

  \author Qi Shao
  $Revision$
  $Date$
*/

template <typename T>
class ParallelInteractionStorage_F : public AParallelInteractionStorage
{

 protected:
  list<CFluidInteraction> m_interactions;
  int m_update_timestamp;

 public:
  ParallelInteractionStorage_F(AParallelParticleArray*);

  virtual ~ParallelInteractionStorage_F(){};

  virtual bool update();
  virtual void calcForces();

  virtual void exchange(){};
  virtual void rebuild(){};
  virtual bool isIn(const vector<int>&){};
  virtual void setTimeStepSize(double dt){};
  virtual void addExIG(AParallelInteractionStorage*){}; // do nothing
  virtual AFieldSlave* generateNewScalarFieldSlave(TML_Comm*,const string&,int,int,int,int){};
  virtual AFieldSlave* generateNewVectorFieldSlave(TML_Comm*,const string&,int,int,int,int){};

  //!< types
  typedef esys::lsm::triplet<int,Vec3,double> ParticleData;

  //!< access functions
  template <typename P>vector<pair<Vec3,P> > forAllInnerInteractionsGetDataWithPos(P (CFluidInteraction::*rdf)() const);
  template <typename P> vector<pair<int,P> > forAllInnerInteractionsGetDataWithID(P (CFluidInteraction::*rdf)() const);
  template <typename P> vector<pair<pair<int,Vec3>,P> > forAllInnerInteractionsGetDataWithIDPos(P (CFluidInteraction::*rdf)() const);
  template <typename P> vector<pair<ParticleData,P> > forAllInnerInteractionsGetDataWithParticle(P (CFluidInteraction::*rdf)() const);
  template <typename P> void forAllInnerInteractionsGet(P& cont,typename P::value_type (CFluidInteraction::*rdf)()const);


};


#include "pis/pi_storage_f.hpp"

#endif //__PARALLEL_INTERACTION_STORAGE_F_H
