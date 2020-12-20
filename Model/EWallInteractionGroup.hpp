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

  k_local=0.0;
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
  }

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
  // warn if trying to apply non-perpendicular force
  double par=F.unit()*this->m_wall->getNormal().unit(); // should be 1 if parallel
  if(par<1.0){
    console.Warning() << "Force not parallel to wall normal in  CEWallInteractionGroup::applyForce : " << F << " vs. " << this->m_wall->getNormal() << "\n";
  }
  // warn if initial value of force is zero, which prevents convergence
  if(F.norm()==0.0){
    console.Warning() << "No force applied to wall in CEWallInteractionGroup::applyForce, which will not converge.  If ramping the force, start with a nonzero value.\n";
  }

  int it=0;
  double d; // distance to move
  double df; // force difference
  double ef; // relative force error (df/|F|)
  Vec3 O_f=F.unit(); // direction of the applied force
  do{
    //std::cerr << "iteration: " << it << std::endl;
    // calculate local stiffness
    k_local=0.0;
    for(typename vector<CElasticWallInteraction<T> >::iterator iter=m_interactions.begin();
        iter!=m_interactions.end();
        iter++){      
      k_local+=iter->getStiffness();
    }
    //std::cerr << "local stiffness: " << k_local <<  std::endl;
    // get global K
    m_k_global=this->m_comm->sum_all(k_local);
    //std::cerr << "global stiffness: " << m_k_global <<  std::endl;
    if(m_k_global>0){
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
      //std::cerr << "local Force: " << F_local <<  std::endl;
      // get global F
      // by component (hack - fix later,i.e. sum_all for Vec3)
      double fgx=this->m_comm->sum_all(F_local.X());
      double fgy=this->m_comm->sum_all(F_local.Y());
      double fgz=this->m_comm->sum_all(F_local.Z());
      Vec3 F_global=Vec3(fgx,fgy,fgz);
      //std::cerr << "global Force: " << F_global <<  std::endl;

      // calc necessary wall movement
      df=(F+F_global)*O_f;
      d=df/m_k_global;
      ef=df/F.norm();
      // move the wall
      this->m_wall->moveBy(d*O_f);
      it++;
    } else {
      d=1e-5;
      ef=1;
      // move the wall
      this->m_wall->moveBy(d*O_f); 
      it++;
    } 
  } while((it<50)&&(ef>1e-3)); // check for convergence
  // warning message if no contact or high error after iteration
  if(ef>1e-3){
    console.Warning() << "applyForce doesn't converge,  global stiffness: " << m_k_global << " applied force: " << F << " wall normal: " << this->m_wall->getNormal() << "\n";
  }
}


template<class T>
ostream& operator<<(ostream& ost,const CEWallInteractionGroup<T>& IG)
{
  ost << "CEWallInteractionGroup" << endl << flush;
  ost << *(IG.m_wall) << endl << flush;
 
  return ost;
}
