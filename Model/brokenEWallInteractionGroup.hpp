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
CEWallInteractionGroup<T>::CEWallInteractionGroup(TML_Comm* comm):AWallInteractionGroup<T>(comm)
{}

/*!
  Constructor for elastic wall interaction group

  \param comm the communicator
  \param wallp a pointer to the wall
  \param param the interaction parameters
*/
template<class T>
CEWallInteractionGroup<T>::CEWallInteractionGroup(TML_Comm* comm,CWall* wallp,const CEWallIGP* I)
  :AWallInteractionGroup<T>(comm)
{
  console.XDebug() << "making CEWallInteractionGroup \n";

  m_k=I->getSpringConst();
  this->m_wall=wallp;
}

template<class T>
void CEWallInteractionGroup<T>::calcForces()
{

  console.XDebug() << "calculating " << m_interactions.size() << " elastic wall forces\n" ;

  for(
    typename vector<CElasticWallInteraction<T> >::iterator it=m_interactions.begin();
    it != m_interactions.end();
    it++
  ){
    it->calcForces();
  }
}

template<class T>
void CEWallInteractionGroup<T>::Update(ParallelParticleArray<T>* PPA)
{

  console.XDebug() << "CEWallInteractionGroup::Update()\n" ;

  console.XDebug()
    << "CEWallInteractionGroup::Update: wall origin = " << this->m_wall->getOrigin()
    << ", wall normal = " << this->m_wall->getNormal() << "\n" ;

  double k_local=0.0;
  // empty particle list first
  m_interactions.erase(m_interactions.begin(),m_interactions.end());
  this->m_inner_count=0;
  // build new particle list
  typename ParallelParticleArray<T>::ParticleListHandle plh=
    PPA->getParticlesAtPlane(this->m_wall->getOrigin(),this->m_wall->getNormal());
  for(typename ParallelParticleArray<T>::ParticleListIterator iter=plh->begin();
      iter!=plh->end();
      iter++){
    bool iflag=PPA->isInInner((*iter)->getPos());
    m_interactions.push_back(CElasticWallInteraction<T>(*iter,this->m_wall,m_k,iflag));
    this->m_inner_count+=(iflag ? 1 : 0);
    if(iflag){
      if(!CParticle::getDo2dCalculations()){
	k_local+=m_k*((*iter)->getRad()); // update local K
      } else {
	k_local+=m_k;
      }
    }
  }
  // get global K
  m_k_global=this->m_comm->sum_all(k_local);

  console.XDebug() << "end CEWallInteractionGroup::Update()\n";
}


/*!
  Apply a given force to the wall. Only forces in the direction of the given force are 
  considered, free movement is assumed in perpendicular directions.

  \param F the force
  \warning Forces not perpendicular to the wall make no sense here, but this is not checked!
*/
template<class T>
void CEWallInteractionGroup<T>::applyForce(const Vec3& F)
{
  int it=0;
  double d;
  Vec3 O_f=F.unit(); // direction of the applied force
  do{
    // calculate local F
    Vec3 F_local=Vec3(0.0,0.0,0.0);
    for (
      typename vector<CElasticWallInteraction<T> >::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++
    ){
      if(iter->isInner()){
        Vec3 f_i=iter->getForce();
        F_local+=(f_i*O_f)*O_f; // add component of f_i in O_f direction
      }
    }
    // get global F
    // by component (hack - fix later,i.e. sum_all for Vec3)
    double fgx=this->m_comm->sum_all(F_local.X());
    double fgy=this->m_comm->sum_all(F_local.Y());
    double fgz=this->m_comm->sum_all(F_local.Z());
    Vec3 F_global=Vec3(fgx,fgy,fgz);

    // calc necessary wall movement
    d=((F+F_global)*O_f)/m_k_global;
    // move the wall
    this->m_wall->moveBy(d*O_f);
    it++;
  } while((it<10)&&(fabs(d)>10e-6)); // check for convergence
}


template<class T>
ostream& operator<<(ostream& ost,const CEWallInteractionGroup<T>& IG)
{
  ost << "CEWallInteractionGroup" << endl << flush;
  ost << *(IG.m_wall) << endl << flush;
 
  return ost;
}
