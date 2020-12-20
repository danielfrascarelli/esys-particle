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

#include "Foundation/console.h"

template<typename P,typename I>
ParallelInteractionStorage_Single<P,I>::ParallelInteractionStorage_Single(AParallelParticleArray* ppa,const typename I::ParameterType& param):TParallelInteractionStorage<I>(ppa)
{
  m_param=param;
}

template<typename P,typename I>
bool ParallelInteractionStorage_Single<P,I>::update()
{
  console.XDebug() << "CDampingGroup<T>::Update()\n" ;
  // empty particle list first
  this->m_interactions.erase(
    this->m_interactions.begin(),
    this->m_interactions.end()
  );
  // build new particle list
  typename ParallelParticleArray<P>::ParticleListHandle plh =
    ((ParallelParticleArray<P>*)this->m_ppa)->getAllParticles();
  for(
    typename ParallelParticleArray<P>::ParticleListIterator iter=plh->begin();
    iter!=plh->end();
    iter++
  ){
    this->m_interactions.push_back(I(*iter,&m_param));
  }
  console.XDebug() << "end CDampingGroup<T>::Update()\n" ;

  return true;
}


template<typename P,typename InteractionType>
void ParallelInteractionStorage_Single<P,InteractionType>::calcForces()
{
  console.Debug()
    << "calculating "
    << this->m_interactions.size()
    << " damping forces\n" ;

  for(
    typename std::list<InteractionType>::iterator it = this->m_interactions.begin();
    it != this->m_interactions.end();
    it++
  ){
    it->calcForces();
  }
}

template<typename P,typename InteractionType>
void
ParallelInteractionStorage_Single<P,InteractionType>::setTimeStepSize(
  double dt
)
{
  console.Debug()
    << "Setting time step size for "
    << this->m_interactions.size()
    << " damping forces\n";

  m_param.setTimeStepSize(dt);
  for (
    typename std::list<InteractionType>::iterator it = this->m_interactions.begin();
    it != this->m_interactions.end();
    it++
  ){
    it->setTimeStepSize(dt);
  }
}
