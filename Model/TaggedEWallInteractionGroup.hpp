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

//----------------------------------------------
//       CEWallInteractionGroup functions
//----------------------------------------------

#include "Foundation/console.h"
#include <iostream>

template<class T>
CTaggedEWallInteractionGroup<T>::CTaggedEWallInteractionGroup(TML_Comm* comm):CEWallInteractionGroup<T>(comm)
{}

/*!
  Constructor for a tagged elastic wall interaction group

  \param comm the communicator
  \param wallp a pointer to the wall
  \param param the interaction parameters
  \param tag the tag of the particles 
  \param mask the mask for the particle tag	
*/
template<class T>
CTaggedEWallInteractionGroup<T>::CTaggedEWallInteractionGroup(TML_Comm* comm,CWall* wallp,const CEWallIGP* I, int tag, int mask)
  :CEWallInteractionGroup<T>(comm, wallp,I)
{
  console.XDebug() << "making CTaggedEWallInteractionGroup \n";

  this->m_tag=tag;
  this->m_mask=mask;		
}

/*!
  Update the interaction group. Checks which particles are close to the wall and have the right
  tag. Except for the tag issue this is copied from the base class (CEWallInteractionGroup)

  \param PPA the array containing the particles
*/
template<class T>
void CTaggedEWallInteractionGroup<T>::Update(ParallelParticleArray<T>* PPA)
{
  console.XDebug() << "CTaggedEWallInteractionGroup::Update()\n" ;

  console.XDebug()
    << "CTaggedEWallInteractionGroup::Update: wall origin = " << this->m_wall->getOrigin()
    << ", wall normal = " << this->m_wall->getNormal() << "\n" ;

  this->k_local=0.0;
  // empty particle list first
  this->m_interactions.erase(this->m_interactions.begin(),this->m_interactions.end());
  this->m_inner_count=0;
  // build new particle list
  typename ParallelParticleArray<T>::ParticleListHandle plh=
    PPA->getParticlesAtPlane(this->m_wall->getOrigin(),this->m_wall->getNormal());
  for(typename ParallelParticleArray<T>::ParticleListIterator iter=plh->begin();
      iter!=plh->end();
      iter++){
	int p_tag=(*iter)->getTag(); // get particle tag
	if ((p_tag & this->m_mask) == (this->m_tag & this->m_mask)){ // check if particles have the right tag
		bool iflag=PPA->isInInner((*iter)->getPos());
		this->m_interactions.push_back(CElasticWallInteraction<T>(*iter,this->m_wall,this->m_k,iflag));
		this->m_inner_count+=(iflag ? 1 : 0);
	}
  }
  console.XDebug() << "found " << this->m_inner_count << " interactions\n";
  
  console.XDebug() << "end CTaggedEWallInteractionGroup::Update()\n";
}

