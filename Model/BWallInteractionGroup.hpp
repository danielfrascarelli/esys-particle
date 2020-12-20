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
//       CBWallInteractionGroup functions
//----------------------------------------------

#include "Foundation/console.h"

template<class T>
CBWallInteractionGroup<T>::CBWallInteractionGroup(TML_Comm* comm):AWallInteractionGroup<T>(comm)
{}

/*!
  Constructor for bonded wall interaction group

  \param comm the communicator
  \param wallp a pointer to the wall
  \param param the interaction parameters
*/
template<class T>
CBWallInteractionGroup<T>::CBWallInteractionGroup(TML_Comm* comm,CWall* wallp,const CBWallIGP* I)
  :AWallInteractionGroup<T>(comm)
{
  console.XDebug() << "making CBWallInteractionGroup \n";

  m_k=I->getSpringConst();
  this->m_wall=wallp;
  m_tag=I->getTag();
  m_mask=I->getMask();
  this->m_inner_count=0;
}

template<class T>
void CBWallInteractionGroup<T>::calcForces()
{

  console.XDebug() << "calculating " << m_bonded_interactions.size() << " bonded wall forces\n" ;
  console.XDebug() << "calculating " << m_elastic_interactions.size() << " elastic wall forces\n" ;

  for (
      typename vector<CBondedWallInteraction<T> >::iterator it=m_bonded_interactions.begin();
      it!=m_bonded_interactions.end();
      it++){
    it->calcForces();
  }
  for(typename vector<CElasticWallInteraction<T> >::iterator it=m_elastic_interactions.begin();
      it!=m_elastic_interactions.end();
      it++){
    it->calcForces();
  }
}


/*!
  Apply a given force to the wall. Only forces in the direction of the given force are 
  considered, free movement is assumed in perpendicular directions.

  \param F the force
*/
template<class T>
void CBWallInteractionGroup<T>::applyForce(const Vec3& F)
{
  console.XDebug() << "CBWallInteractionGroup<T>::applyForce: F = " << F << "\n";
  // calculate local K
  
  double K=0.0;
  for (
      typename vector<CBondedWallInteraction<T> >::iterator it=m_bonded_interactions.begin();
      it!=m_bonded_interactions.end();
      it++){
    if(it->isInner()){
      K+=it->getStiffness();
    }
  }
  for(typename vector<CElasticWallInteraction<T> >::iterator it=m_elastic_interactions.begin();
      it!=m_elastic_interactions.end();
      it++){
    if(it->isInner()){
      K+=it->getStiffness();
    }
  }
  // get global K
  double K_global=this->m_comm->sum_all(K);

  int it=0;
  double d;
  Vec3 O_f=F.unit(); // direction of the applied force
  console.XDebug() << "CBWallInteractionGroup<T>::applyForce: unitF = " << O_f << "\n";
  do{
    // calculate local F
    Vec3 F_local=Vec3(0.0,0.0,0.0);
    // bonded itneractions
    for(
      typename vector<CBondedWallInteraction<T> >::iterator iter=m_bonded_interactions.begin();
      iter!=m_bonded_interactions.end();
      iter++
    ){
      if(iter->isInner()){
        Vec3 f_i=iter->getForce();
        F_local+=(f_i*O_f)*O_f; // add component of f_i in O_f direction
      }
    }
    // elastic interactions
    for(
      typename vector<CElasticWallInteraction<T> >::iterator iter=m_elastic_interactions.begin();
      iter!=m_elastic_interactions.end();
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
    console.XDebug()
      << "CBWallInteractionGroup<T>::applyForce: iteration " << it << ", d = " << fabs(d) << "\n";
    
    // move the wall
    console.XDebug()
      << "CBWallInteractionGroup<T>::applyForce: moving wall by " << d*O_f << "\n";    
    this->m_wall->moveBy(d*O_f);
    it++;
  } while((it<10)&&(fabs(d)>10e-6)); // check for convergence
  console.XDebug()
    << "CBWallInteractionGroup<T>::applyForce: d = " << fabs(d)
    << ", num iterations = " << it << "\n";
}

/*!
  Update interactions from an existing parallel particle array

  \param PPA a pointer to the particle array
*/
template<class T>
void CBWallInteractionGroup<T>::Update(ParallelParticleArray<T>* PPA)
{

  console.XDebug() << "CBWallInteractionGroup::Update()\n" ;

  // -- bonded interactions --
  // empty particle list first
  m_bonded_interactions.erase(m_bonded_interactions.begin(),m_bonded_interactions.end());
  this->m_inner_count=0;
  // build new particle list
  typename ParallelParticleArray<T>::ParticleListHandle plh=
    PPA->getParticlesAtPlane(this->m_wall->getOrigin(),this->m_wall->getNormal());
  for(typename ParallelParticleArray<T>::ParticleListIterator iter=plh->begin();
      iter!=plh->end();
      iter++){
    if(((*iter)->getTag() & m_mask )== (m_tag & m_mask)){ // if tagged->bonded ,add to bonded interaction
      bool iflag=PPA->isInInner((*iter)->getPos());
      m_bonded_interactions.push_back(CBondedWallInteraction<T>(*iter,this->m_wall,m_k,iflag));
      this->m_inner_count+=(iflag ? 1 : 0);
    }
  }
  // -- elastic interactions --
  // empty particle list first
  m_elastic_interactions.erase(m_elastic_interactions.begin(),m_elastic_interactions.end());
  // build new particle list
  for(typename ParallelParticleArray<T>::ParticleListIterator iter=plh->begin();
      iter!=plh->end();
      iter++){
    if(((*iter)->getTag() & m_mask )!=(m_tag & m_mask )){ // only add to elastic if not bonded/tagged
      bool iflag=PPA->isInInner((*iter)->getPos());
      m_elastic_interactions.push_back(CElasticWallInteraction<T>(*iter,this->m_wall,m_k,iflag));
      this->m_inner_count+=(iflag ? 1 : 0);
    }
  }
  console.XDebug() << "end CBWallInteractionGroup::Update()\n";
}

template<class T>
ostream& operator<<(ostream& ost,const CBWallInteractionGroup<T>& IG)
{
  ost << "CBWallInteractionGroup" << endl << flush;
  ost << *(IG.m_wall) << endl << flush;
 
  return ost;
}
