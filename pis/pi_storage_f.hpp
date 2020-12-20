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


#include "Foundation/vec3.h"
#include "ppa/src/pp_array.h"
#include "Model/FluidInteraction.h"
#include "Fields/ScalarFluidInteractionFieldSlave.h"
#include "Fields/VectorFluidInteractionFieldSlave.h"

typedef esys::lsm::triplet<int,Vec3,double> ParticleData;



template<typename T>
ParallelInteractionStorage_F<T>::ParallelInteractionStorage_F(AParallelParticleArray* ppa):AParallelInteractionStorage(ppa)
{
  m_update_timestamp=0;
}


/*!
  Update interactions. Do full dynamic search.
*/
template <typename T>
bool ParallelInteractionStorage_F<T>::update()
{
  console.XDebug() << "ParallelInteractionStorage_F::Update\n";
  int count_l=0;

  // clean out old interactions
  this->m_interactions.clear();

  // get list  of pairs from m_ppa
  typename ParallelParticleArray<T>::ParticleCellListHandle plh =
    ((ParallelParticleArray<T>*)this->m_ppa)->getParticleCellPairList();

  // generate interactions from pairs
  for(typename ParallelParticleArray<T>::ParticleCellListIterator iter=plh->begin();
      iter!=plh->end();
      iter++){
    this->m_interactions.push_back(CFluidInteraction(iter->first,iter->second));
    count_l++;
  }

  console.XDebug() << "added " << count_l << " pairs to FluidInteractions\n";
  console.XDebug() << "end ParallelInteractionStorage_F::Update\n";
  return true;
}


template <typename T>
void ParallelInteractionStorage_F<T>::calcForces()
{
  console.Debug()
    << "calculating "
    << this->m_interactions.size()
    << " fluid forces\n" ;
  int n=0;
  for(
    typename std::list<CFluidInteraction>::iterator it = this->m_interactions.begin();
    it != this->m_interactions.end();
    it++
  ){
    it->calcForces();
    n++;
  }
}


/*!
  For all interactions with the particle in the inner area of the ntable call a
  function reading a value and return the results in a vector of <p_pos,vaule> pairs

  \param rdf the function
*/
template <typename T>
template <typename P>
vector<pair<Vec3,P> > ParallelInteractionStorage_F<T>::forAllInnerInteractionsGetDataWithPos(P (CFluidInteraction::*rdf)() const)
{
  vector<pair<Vec3,P> > res;

  for(list<CFluidInteraction>::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++){
    Vec3 pos=iter->getParticlePos();
    if(m_ppa->isInInner(pos)) res.push_back(make_pair(pos,((*iter).*rdf)()));
  }
  return res;
}

/*!
  For all interactions with the particle in the inner area of the ntable call a
  function reading a value and return the results in a vector of <p_id,vaule> pairs

  \param rdf the function
*/
template <typename T>
template <typename P>
vector<pair<int,P> > ParallelInteractionStorage_F<T>::forAllInnerInteractionsGetDataWithID(P (CFluidInteraction::*rdf)() const)
{
  vector<pair<int,P> > res;

  for(list<CFluidInteraction>::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++){
    int id=iter->getParticleID();
    Vec3 pos=iter->getParticlePos();
    if(m_ppa->isInInner(pos)) res.push_back(make_pair(id,((*iter).*rdf)()));
  }
  return res;
}

/*!
  For all interactions with the particle in the inner area of the ntable call a
  function reading a value and return the results in a vector of <<p_id,p_pos>,value> groups

  \param rdf the function
*/
template <typename T>
template <typename P>
vector<pair<pair<int,Vec3>,P> >
ParallelInteractionStorage_F<T>::forAllInnerInteractionsGetDataWithIDPos(P (CFluidInteraction::*rdf)() const)
{
  vector<pair<pair<int,Vec3>,P> > res;

  for(list<CFluidInteraction>::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++){
    if(m_ppa->isInInner(iter->getParticlePos())) {
      int id=iter->getParticleID();
      Vec3 pos=iter->getParticlePos();
      res.push_back(pair<pair<int,Vec3>,P>(pair<int,Vec3>(id,pos),((*iter).*rdf)()));
    }
  }

  return res;
}


/*!
  For all interactions with the particle in the inner area of the ntable call a
  function reading a value and return the results in a vector of <<p_id,p_pos,p_r>,value> groups

  \param rdf the function
*/
template <typename T>
template <typename P>
vector<pair<ParticleData,P> >
ParallelInteractionStorage_F<T>::forAllInnerInteractionsGetDataWithParticle(P (CFluidInteraction::*rdf)() const)
{
  vector<pair<ParticleData,P> > res;

  for(list<CFluidInteraction>::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++){
    if(m_ppa->isInInner(iter->getParticlePos())) {
      const ParticleData data = iter->getParticleData();
      res.push_back(pair<ParticleData,P>(data,((*iter).*rdf)()));
    }
  }

  return res;
}


/*!
  For all interactions with the particle in the inner area of the ntable call a
  function reading a value and return the results in a container
  particle ids

  \cont the container
  \param rdf the function
*/
template <typename T>
template <typename P>
void ParallelInteractionStorage_F<T>::forAllInnerInteractionsGet(P& cont,typename P::value_type (CFluidInteraction::*rdf)()const)
{
  for(typename list<CFluidInteraction>::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++){
    Vec3 pos=iter->getParticlePos();
    if(m_ppa->isInInner(pos)) cont.push_back(((*iter).*rdf)());
  }
}


