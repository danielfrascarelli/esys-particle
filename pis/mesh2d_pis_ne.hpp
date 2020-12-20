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

/*!
  constructor

  \param mesh_p pointer to the 2d mesh
  \param ppa_p pointer to the particle array
  \param param 
*/
template<class ParticleType,class IType>
Mesh2D_PIS_NE<ParticleType,IType>::Mesh2D_PIS_NE(Mesh2D* mesh_p,ParallelParticleArray<ParticleType>* ppa_p,typename IType::ParameterType param)
  :Mesh2D_PIS<ParticleType>(mesh_p,ppa_p)
{
  m_param=param;
  this->m_update_timestamp=0;
}

/*!
  destructor
*/
template<class ParticleType,class IType>
Mesh2D_PIS_NE<ParticleType,IType>::~Mesh2D_PIS_NE()
{}

/*!
  Check if an interaction is in this PIS. The first 2 values in the vector
  are expected to be the edge/corner (v[0]) and particle (v[1]) ids, the 3rd 
  an indicator if edge (v[2]==0)or corner (v[2]==1) interaction. 
  If there is no 3rd value or it is not in [0,1], "false" is returned.

  \param v vector of particle ids
  \warning log(N)
*/
template <class ParticleType,class IType> 
bool Mesh2D_PIS_NE<ParticleType,IType>::isIn(const std::vector<int>& v)
{
  bool res=false;
  
  if(v.size()<3){
    res=false;
  } else {
    switch (v[2]){
    case 0: res=m_edge_int_set.find(make_pair(v[0],v[1]))!=m_edge_int_set.end(); break;
    case 1: res=m_corner_int_set.find(make_pair(v[0],v[1]))!=m_corner_int_set.end(); break;
    default: console.Error() << "wrong value in argument of Mesh2D_PIS::isIn !!\n"; break;
   }
  }

  return res;
}

/*!
  calculate all the forces
*/
template<class ParticleType,class IType>
void Mesh2D_PIS_NE<ParticleType,IType>::calcForces()
{
  console.XDebug() << "Mesh2D_PIS_NE calculating " << m_edge_interactions.size() << " line forces , " 
		   << m_corner_interactions.size() << "corner forces\n";

  // calculate forces for edge interactions
  for(typename std::vector<typename IType::EdgeIntType>::iterator tri_iter=m_edge_interactions.begin();
      tri_iter!=m_edge_interactions.end();
      tri_iter++){
    tri_iter->calcForces();
  }
  // calculate forces for corner interactions 
  for(typename std::vector<typename IType::CornerIntType>::iterator corner_iter=m_corner_interactions.begin();
      corner_iter!=m_corner_interactions.end();
      corner_iter++){
    corner_iter->calcForces();
  }
}

/*!
  update the interactions 
*/
template<class ParticleType,class IType>
bool Mesh2D_PIS_NE<ParticleType,IType>::update()
{
  console.XDebug() << "Mesh2D_PIS_NE::update\n";
  bool res=false;
  //int count_edge=0;
  //int count_tri=0;

  if(this->m_update_timestamp != this->m_ppa->getTimeStamp()){// m_ppa rebuild since last update 
    console.XDebug() << "Mesh2D_PIS_NE doing update\n";
    // clean out old interactions
    m_edge_interactions.clear();
    m_corner_interactions.clear();
    m_edge_int_set.clear();
    m_corner_int_set.clear();
    // -- get edge interactions
    // for all edges
    for(
      Mesh2D::edge_iterator ed_iter = this->m_mesh->edges_begin();
      ed_iter != this->m_mesh->edges_end();
      ed_iter++
    ){
      typename ParallelParticleArray<ParticleType>::ParticleListHandle plh =
        ((ParallelParticleArray<ParticleType>*)(this->m_ppa))->getParticlesNearEdge(&(*ed_iter));
      for (
        typename ParallelParticleArray<ParticleType>::ParticleListIterator p_iter=plh->begin();
        p_iter!=plh->end();
        p_iter++
      ){
        bool iflag = this->m_ppa->isInInner((*p_iter)->getPos());
        m_edge_interactions.push_back(typename IType::EdgeIntType(*p_iter,&(*ed_iter),m_param,iflag));
        //m_particle_id_set.insert((*p_iter)->getID());
      }
    }
    // --- get corner interactions
    for (
      Mesh2D::corner_iterator co_iter = this->m_mesh->corners_begin();
      co_iter != this->m_mesh->corners_end();
      co_iter++
    ){
      typename ParallelParticleArray<ParticleType>::ParticleListHandle plh=
        ((ParallelParticleArray<ParticleType>*)(this->m_ppa))->getParticlesNearPoint(co_iter->getPos());
      for (
        typename ParallelParticleArray<ParticleType>::ParticleListIterator p_iter=plh->begin();
        p_iter!=plh->end();
        p_iter++
      ){
        bool iflag = this->m_ppa->isInInner((*p_iter)->getPos());
        m_corner_interactions.push_back(typename IType::CornerIntType(*p_iter,&(*co_iter),m_param,iflag));
        //m_particle_id_set.insert((*p_iter)->getID());
      }
    }
    // set timestamp
    this->m_update_timestamp = this->m_ppa->getTimeStamp();
  }
  console.XDebug() << "end  ElasticMesh2DIG<T>::Update\n";

  return res;
}

