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
//       CSoftBWallInteractionGroup functions
//----------------------------------------------

#include "Foundation/console.h"

template<class T>
CSoftBWallInteractionGroup<T>::CSoftBWallInteractionGroup(TML_Comm* comm):AWallInteractionGroup<T>(comm)
{}

/*!
  Constructor for bonded wall interaction group with direction dependend elasticity

  \param comm the communicator
  \param wallp a pointer to the wall
  \param param the interaction parameters
*/
template<class T>
CSoftBWallInteractionGroup<T>::CSoftBWallInteractionGroup(TML_Comm* comm,CWall* wallp, const CSoftBWallIGP* I)
  : AWallInteractionGroup<T>(comm)
{
  console.XDebug() << "making CSoftBWallInteractionGroup \n";

  m_normalK=I->getNormalK();
  m_shearK=I->getShearK();
  this->m_wall=wallp;
  m_tag=I->getTag();
  m_mask=I->getMask();
  m_scaling=I->getScaling();
//  console.XDebug() << "kx, ky, kz: " << m_kx << ","<< m_ky << ","<< m_kz << "\n";
}

template<class T>
void CSoftBWallInteractionGroup<T>::calcForces()
{

  console.XDebug() << "calculating " << m_interactions.size() << " soft bonded wall forces\n" ;

  for(
    typename vector<CSoftBondedWallInteraction<T> >::iterator it=m_interactions.begin();
    it != m_interactions.end();
    it++
  ){
    it->calcForces();
  }
}

/*!
  Apply a given force to the wall. Only forces in the direction of the given force are 
  considered, free movement is assumed in perpendicular directions.

  \param F the force
*/
template<class T>
void CSoftBWallInteractionGroup<T>::applyForce(const Vec3& F)
{
  console.XDebug() << "CSoftBWallInteractionGroup<T>::applyForce: F = " << F << "\n";
  // calculate local K

  double K=0.0;
  for (
      typename vector<CSoftBondedWallInteraction<T> >::iterator it=m_interactions.begin();
      it!=m_interactions.end();
      it++){
    if(it->isInner()){
      K+=it->getStiffness();
      console.XDebug() << "CSoftBWallInteractionGroup<T>::applyForce: K = " << K << "\n";
    }
  }
  // get global K
  double K_global=this->m_comm->sum_all(K);
  console.XDebug() << "CSoftBWallInteractionGroup<T>::applyForce: K_global = " << K_global << "\n";

  int it=0;
  double d;
  Vec3 O_f=F.unit(); // direction of the applied force
  console.XDebug() << "CSoftBWallInteractionGroup<T>::applyForce: unitF = " << O_f << "\n";
  do{
    // calculate local F
    Vec3 F_local=Vec3(0.0,0.0,0.0);
    // bonded itneractions
    for(
      typename vector<CSoftBondedWallInteraction<T> >::iterator iter=m_interactions.begin();
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


template<class T>
void CSoftBWallInteractionGroup<T>::Update(ParallelParticleArray<T>* PPA)
{

  console.XDebug() << "CSoftBWallInteractionGroup::Update()\n" ;

  // empty particle list first
  m_interactions.erase(m_interactions.begin(),m_interactions.end());
  // build new particle list
  typename ParallelParticleArray<T>::ParticleListHandle plh=
    PPA->getParticlesAtPlane(this->m_wall->getOrigin(),this->m_wall->getNormal());
  for(typename ParallelParticleArray<T>::ParticleListIterator iter=plh->begin();
      iter!=plh->end();
      iter++){
    if(((*iter)->getTag() & m_mask) ==(m_tag & m_mask)){
      bool iflag=PPA->isInInner((*iter)->getPos());
      m_interactions.push_back(CSoftBondedWallInteraction<T>(*iter,this->m_wall,m_normalK,m_shearK,m_scaling,iflag));
    }
  }
  console.XDebug() << "end CSoftBWallInteractionGroup::Update()\n";
}

template<class T>
ostream& operator<<(ostream& ost,const CSoftBWallInteractionGroup<T>& IG)
{
  ost << "CBWallInteractionGroup" << endl << flush;
  ost << *(IG.m_wall) << endl << flush;
 
  return ost;
}
