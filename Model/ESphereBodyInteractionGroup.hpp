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
//       CESphereBodyInteractionGroup functions
//----------------------------------------------

#include "Foundation/console.h"
#include <iostream>

template<class T>
CESphereBodyInteractionGroup<T>::CESphereBodyInteractionGroup(TML_Comm* comm):ASphereBodyInteractionGroup<T>(comm)
{}

/*!
  Constructor for elastic sphere body interaction group

  \param comm the communicator
  \param spherep a pointer to the sphere body
  \param param the interaction parameters
*/
template<class T>
CESphereBodyInteractionGroup<T>::CESphereBodyInteractionGroup(TML_Comm* comm,CSphereBody* spherep,const CESphereBodyIGP* I)
  :ASphereBodyInteractionGroup<T>(comm)
{
  console.XDebug() << "making CESphereBodyInteractionGroup \n";

  m_k=I->getSpringConst();
  this->m_sphere=spherep;
}

template<class T>
void CESphereBodyInteractionGroup<T>::calcForces()
{

  console.XDebug() << "calculating " << m_interactions.size() << " elastic sphere body forces\n" ;

  for(
    typename vector<CElasticSphereBodyInteraction<T> >::iterator it=m_interactions.begin();
    it != m_interactions.end();
    it++
  ){
    it->calcForces();
  }
}

template<class T>
void CESphereBodyInteractionGroup<T>::Update(ParallelParticleArray<T>* PPA)
{

  console.XDebug() << "CESphereBodyInteractionGroup::Update()\n" ;

  console.XDebug()
    << "CESphereBodyInteractionGroup::Update: sphere body origin = " << this->m_sphere->getCentre()
    << ", sphere body radius = " << this->m_sphere->getRadius() << "\n" ;

  k_local=0.0;
  // empty particle list first
  m_interactions.erase(m_interactions.begin(),m_interactions.end());
  this->m_inner_count=0;
  // build new particle list
  typename ParallelParticleArray<T>::ParticleListHandle plh=
    PPA->getParticlesNearSphere(this->m_sphere->getCentre(),this->m_sphere->getRadius());
  for(typename ParallelParticleArray<T>::ParticleListIterator iter=plh->begin();
      iter!=plh->end();
      iter++){
    bool iflag=PPA->isInInner((*iter)->getPos());
    m_interactions.push_back(CElasticSphereBodyInteraction<T>(*iter,this->m_sphere,m_k,iflag));
    this->m_inner_count+=(iflag ? 1 : 0);
  }

  console.XDebug() << "end CESphereBodyInteractionGroup::Update()\n";
}


/*!
  Apply a given force to the sphere body. Only forces in the direction of the given force are 
  considered, free movement is assumed in perpendicular directions.

  \param F the force
*/
template<class T>
void CESphereBodyInteractionGroup<T>::applyForce(const Vec3& F)
{
  // warn if initial value of force is zero, which prevents convergence
  if(F.norm()==0.0){
    console.Warning() << "No force applied to sphere body in CESphereBodyInteractionGroup::applyForce, which will not converge.  If ramping the force, start with a nonzero value.\n";
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
    for(typename vector<CElasticSphereBodyInteraction<T> >::iterator iter=m_interactions.begin();
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
          typename vector<CElasticSphereBodyInteraction<T> >::iterator iter=m_interactions.begin();
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

      // calc necessary sphere movement
      df=(F+F_global)*O_f;
      d=df/m_k_global;
      ef=df/F.norm();
      // move the sphere body
      this->m_sphere->moveBy(d*O_f);
      it++;
    } else {
      d=1e-5;
      ef=1;
      // move the sphere body
      this->m_sphere->moveBy(d*O_f); 
      it++;
    } 
  } while((it<50)&&(ef>1e-3)); // check for convergence
  // warning message if no contact or high error after iteration
  if(ef>1e-3){
    console.Warning() << "applyForce doesn't converge,  global stiffness: " << m_k_global << " applied force: " << F << "\n";
  }
}


template<class T>
ostream& operator<<(ostream& ost,const CESphereBodyInteractionGroup<T>& IG)
{
  ost << "CESphereBodyInteractionGroup" << endl << flush;
  ost << *(IG.m_sphere) << endl << flush;
 
  return ost;
}
