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
//       CViscWallIG functions
//----------------------------------------------

#include "Foundation/console.h"

template<class T>
CViscWallIG<T>::CViscWallIG(TML_Comm* comm):AWallInteractionGroup<T>(comm)
{
}

/*!
  Constructor for wall interaction group with viscous drag

  \param comm the communicator
  \param wallp a pointer to the wall
  \param param the interaction parameters
*/
template<class T>
CViscWallIG<T>::CViscWallIG(TML_Comm* comm,CWall* wallp,const CVWallIGP* I)
  :AWallInteractionGroup<T>(comm)
{
  console.XDebug() << "making CViscWallIG pos \n";

  m_k=I->getSpringConst();
  this->m_wall=wallp;
  m_tag=I->getTag();
  m_nu=I->getNu();
  this->m_inner_count=0;
}

template<class T>
void CViscWallIG<T>::calcForces()
{

  console.XDebug() << "calculating " << m_visc_interactions.size() << " viscous wall forces\n" ;
  console.XDebug() << "calculating " << m_elastic_interactions.size() << " elastic wall forces\n" ;

  for(
    typename vector<CViscWallInteraction<T> >::iterator it=m_visc_interactions.begin();
    it!=m_visc_interactions.end();
    it++
  ){
    it->calcForces();
  }
  for(
    typename vector<CElasticWallInteraction<T> >::iterator it=m_elastic_interactions.begin();
    it!=m_elastic_interactions.end();
    it++
  ){
    it->calcForces();
  }
}

/*!
  Set velocity of the wall. Only sets m_vel of the wall, doesn't affect position updates.

  \param V the velocity
*/
template<class T>
void CViscWallIG<T>::setVelocity(const Vec3& V)
{
  this->m_wall->setVel(V);
}

/*!
  Apply a given force to the wall. Only forces in the direction of the given force are 
  considered, free movement is assumed in perpendicular directions.

  \param F the force
*/
template<class T>
void CViscWallIG<T>::applyForce(const Vec3& F)
{
  // calculate local K
  double K=this->m_inner_count*m_k;
  // get global K
  double K_global=this->m_comm->sum_all(K);

  int it=0;
  double d;
  Vec3 O_f=F.unit(); // direction of the applied force
  do{
    // calculate local F
    Vec3 F_local=Vec3(0.0,0.0,0.0);
    // viscous interactions
    for (
      typename vector<CViscWallInteraction<T> >::iterator iter=m_visc_interactions.begin();
        iter != m_visc_interactions.end();
        iter++
      ){
      if(iter->isInner()){
        Vec3 f_i=iter->getForce();
        F_local+=(f_i*O_f)*O_f; // add component of f_i in O_f direction
      }
    }
    // elastic interactions
    for (
      typename vector<CElasticWallInteraction<T> >::iterator iter=m_elastic_interactions.begin();
      iter != m_elastic_interactions.end();
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
    d=((F+F_global)*O_f)/K_global;
    // move the wall
    this->m_wall->moveBy(d*O_f);
    it++;
  } while((it<10)&&(fabs(d)>10e-6)); // check for convergence
}

/*!
  Update interactions from an existing parallel particle array

  \param PPA a pointer to the particle array
*/
template<class T>
void CViscWallIG<T>::Update(ParallelParticleArray<T>* PPA)
{

  console.XDebug() << "CViscWallIG::Update()\n" ;

  // -- bonded interactions --
  // empty particle list first
  m_visc_interactions.erase(m_visc_interactions.begin(),m_visc_interactions.end());
  this->m_inner_count=0;
  // build new particle list
  typename ParallelParticleArray<T>::ParticleListHandle plh=
    PPA->getParticlesAtPlane(this->m_wall->getOrigin(),this->m_wall->getNormal());
  for(typename ParallelParticleArray<T>::ParticleListIterator iter=plh->begin();
      iter!=plh->end();
      iter++){
    if((*iter)->getTag()==m_tag){// if tagged -> apply viscous drag, i.e. add to viscous interaction
      bool iflag=PPA->isInInner((*iter)->getPos());
      m_visc_interactions.push_back(CViscWallInteraction<T>(*iter,this->m_wall,m_nu,iflag));
      this->m_inner_count+=(iflag ? 1 : 0);
    }
  }
  // -- elastic interactions --
  // empty particle list first
  m_elastic_interactions.erase(m_elastic_interactions.begin(),m_elastic_interactions.end());
  // build new particle list
  for(typename ParallelParticleArray<T>::ParticleListIterator iter=plh->begin();
      iter!=plh->end();
      iter++){ // always add to elastic group, even if in viscous
    bool iflag=PPA->isInInner((*iter)->getPos());
    m_elastic_interactions.push_back(CElasticWallInteraction<T>(*iter,this->m_wall,m_k,iflag));
    this->m_inner_count+=(iflag ? 1 : 0);
  }
  console.XDebug() << "end CViscWallIG::Update()\n";
}

template<class T>
ostream& operator<<(ostream& ost,const CViscWallIG<T>& IG)
{
  ost << "CViscWallIG" << endl << flush;
  ost << *(IG.m_wall) << endl << flush;
 
  return ost;
}
