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

template<class ParticleType,class IType>
const int TriMesh_PIS_EB<ParticleType,IType>::m_exchg_tag=44;

/*!
  constructor

  \param mesh_p a pointer to the triangle mesh
  \param ppa_p a pointer to the particle array
  \param param the interaction parameters
*/
template<class ParticleType,class IType>
TriMesh_PIS_EB<ParticleType,IType>::TriMesh_PIS_EB(TriMesh* mesh_p,ParallelParticleArray<ParticleType>* ppa_p,typename IType::ParameterType param)
  : TriMesh_PIS<ParticleType>(mesh_p,ppa_p),m_comm(ppa_p->getComm())
{
  console.XDebug() << "TriMesh_PIS_EB constructor\n";
  m_param=param;
  this->m_update_timestamp = 0;
}
   
template <class ParticleType,class IType> 
void TriMesh_PIS_EB<ParticleType,IType>::setTimeStepSize(double dt)
{
}

/*!
  Check if an interaction is in this PIS. The first 2 values in the vector
  are expected to be the tri/edge/corner (v[0]) and particle (v[1]) ids, the 3rd 
  an indicator if tri (v[2]==0),edge (v[2]==1)or corner (v[2]==2) interaction. 
  If there is no 3rd value or it is not 0 (tri), "false" is returned.

  \param v vector of particle ids
  \warning log(N)
*/
template <class ParticleType,class IType> 
bool TriMesh_PIS_EB<ParticleType,IType>::isIn(const std::vector<int>& v)
{
  bool res=false;
  
  if(v.size()<3){
    res=false;
  } else {
    switch (v[2]){
    case 0: res=m_tri_int_set.find(make_pair(v[0],v[1]))!=m_tri_int_set.end(); break;
    default: console.Error() << "wrong value in argument of TriMesh_PIS::isIn !!\n"; break;
   }
  }

  return res;
}

/*!
  calculate all the forces
*/
template<class ParticleType,class IType>
void TriMesh_PIS_EB<ParticleType,IType>::calcForces()
{
  console.XDebug() << "TriMesh_PIS_EB calculating " << m_triangle_interactions.size() << " triangle forces\n";

  // calculate forces for triangle interactions
  for(typename list<typename IType::TriIntType>::iterator tri_iter=m_triangle_interactions.begin();
      tri_iter!=m_triangle_interactions.end();
      tri_iter++){
    tri_iter->calcForces();
  }
}

/*!
 */
template<class ParticleType,class IType>
bool TriMesh_PIS_EB<ParticleType,IType>::update()
{
  console.XDebug() << "TriMesh_PIS_EB::update on node " << m_comm.rank() << "\n"; 
  bool res=false;

  typename list<typename IType::TriIntType>::iterator iter=m_triangle_interactions.begin();
  while(iter!=m_triangle_interactions.end()){
    if(iter->broken()){
      res=true;
      typename list<typename IType::TriIntType>::iterator er_iter=iter;
      iter++;
      // remove ids from map
      m_tri_int_set.erase(make_pair(er_iter->getTid(),er_iter->getPid()));
      // remove interaction
      m_triangle_interactions.erase(er_iter);
    } else {
      iter++;
    }
  }

  console.XDebug() << "end TriMesh_PIS_EB::update on node " << m_comm.rank() << "\n"; 
  return res;
}
   
/*!
  helper function to do the actual shifting of values in exchange()

  \param dim dimension, 0->x, 1->y, 2->z
  \param dir direction, 1->up, -1->down
*/
template<class ParticleType,class IType>
void TriMesh_PIS_EB<ParticleType,IType>::exchange_boundary(int dim,int dir)
{
  console.XDebug() << "TriMesh_PIS_EB::exchange_boundary(" << dim << "," << dir << ") at node " << m_comm.rank() << "\n";
  
  std::set<int> bdry_ids;
  std::vector<typename IType::TriIntType> recv_tri_buffer;
  std::vector<typename IType::TriIntType> send_tri_buffer;

  // --- setup data to send ---
  // get boundary
  bdry_ids = this->m_ppa->getBoundarySlabIds(dim,dir);
  // for all interactions
  for(typename list<typename IType::TriIntType>::iterator iter=m_triangle_interactions.begin();
      iter!=m_triangle_interactions.end();
      iter++){
    int pid=iter->getPid(); 
    if(bdry_ids.find(pid)!=bdry_ids.end()) { // if particle in interaction is in bdry -> put in buffer
      send_tri_buffer.push_back(*iter);
    }  
  }
  // ---- shift ----
  m_comm.shift_cont_packed(send_tri_buffer,recv_tri_buffer,dim,dir,m_exchg_tag);
  // insert received data
  for(typename std::vector<typename IType::TriIntType>::iterator iter=recv_tri_buffer.begin();
      iter!=recv_tri_buffer.end();
      iter++){
    tryInsert(*iter);
  }
  
  console.XDebug() << "end TriMesh_PIS_EB::exchange_boundary\n";
}

/*!
 */
template<class ParticleType,class IType>
void TriMesh_PIS_EB<ParticleType,IType>::exchange()
{
  console.XDebug() << "TriMesh_PIS_EB::exchange\n";
  for(int i=0;i<3;i++){
    if(m_comm.get_dim(i)>1){
      // -- up --
      exchange_boundary(i,1);
      // -- down --
      exchange_boundary(i,-1);
    }
  }
  console.XDebug() << "end TriMesh_PIS_EB::exchange\n";
}
   
/*!
  Rebuild interactions after moving particles or interactions between processes. Set 
  particle pointers accordig to particle IDs and remove interactionw which include
  unavailable particles.
*/
template<class ParticleType,class IType>
void TriMesh_PIS_EB<ParticleType,IType>::rebuild()
{
  console.XDebug() << "TriMesh_PIS_EB::rebuild on node " << m_comm.rank() << "\n"; 
  ParallelParticleArray<ParticleType>* t_ppa =
    (ParallelParticleArray<ParticleType>*)(this->m_ppa); // should be a dynamic_cast
  // for all triangle interactions
  typename list<typename IType::TriIntType>::iterator ti_iter=m_triangle_interactions.begin();
  while(ti_iter!=m_triangle_interactions.end()){
    int pid=ti_iter->getPid();
    ParticleType *part_p=t_ppa->getParticlePtrByIndex(pid);
    if(part_p!=NULL) { // particle is available in m_ppa -> setup pointer in interaction
      ti_iter->setPP(part_p);
      Triangle *tri_p = this->m_mesh->getTriangleById(ti_iter->getTid());
      ti_iter->setTP(tri_p);
      ti_iter++;
    } else { // particle is not available in m_ppa -> remove interaction
      const typename list<typename IType::TriIntType>::iterator er_iter=ti_iter;
      ti_iter++;
      m_tri_int_set.erase(make_pair(er_iter->getTid(),pid));
      m_triangle_interactions.erase(er_iter); 
    }
  }

  console.XDebug() << "end TriMesh_PIS_EB::rebuild on node " << m_comm.rank() << "\n"; 
}
   
/*!
  Insert interactions newly created from particle Ids and parameters. If insertion is impossible 
  because the interaction is already in, or one of the particles is not in the associated PPA 
  nothing happens.
  Check if an interaction is in this PIS. The first 2 values in the vector
  are expected to be the tri/edge/corner (pids[0]) and particle (pids[1]) ids, the 3rd 
  an indicator if tri (pids[2]==0),edge (pids[2]==1)or corner (pids[2]==2) interaction. 
  If there is no 3rd value or it is not in [0,1,2], nothing happens.

  \param pids the particle Ids
*/
template<class ParticleType,class IType>
void TriMesh_PIS_EB<ParticleType,IType>::tryInsert(const std::vector<int>& pids)
{
  ParallelParticleArray<ParticleType>* t_ppa =
    (ParallelParticleArray<ParticleType>*)(this->m_ppa); // should be dynamic_cast

  if(pids.size()<3){
    bool is_in=isIn(pids); // interaction already in
    ParticleType *part_p=t_ppa->getParticlePtrByIndex(pids[1]);
    if((!is_in) && (part_p!=NULL)){ 
      switch (pids[2]){
        case 0: {
          Triangle *tri_p = this->m_mesh->getTriangleById(pids[0]);
          if (tri_p!=NULL){
            m_triangle_interactions.push_back(
              typename IType::TriIntType(
                part_p,
                tri_p,
                m_param,
                this->m_ppa->isInInner(part_p->getPos())
              )
            );
            m_tri_int_set.insert(make_pair(pids[0],pids[1]));
          }
        } break;
        default : {
          console.Error()
            << "Error: pids[2]= " << pids[2]
            << "\n"; // Can't happen !! 
        }
      }
    }
  }
}

/*!
 */
template<class ParticleType,class IType>
void TriMesh_PIS_EB<ParticleType,IType>::tryInsert(const typename IType::TriIntType& In)
{
  ParallelParticleArray<ParticleType>* t_ppa =
    (ParallelParticleArray<ParticleType>*)(this->m_ppa);
  // check if interaction is already in 
  bool is_in=(m_tri_int_set.find(make_pair(In.getTid(),In.getPid()))!=m_tri_int_set.end());
  ParticleType *part_p=t_ppa->getParticlePtrByIndex(In.getPid());
  if((!is_in) && (part_p!=NULL)){
    m_triangle_interactions.push_back(In);
    m_tri_int_set.insert(make_pair(In.getTid(),In.getPid()));
  }
}

/*!
 */
template<class ParticleType,class IType>
void TriMesh_PIS_EB<ParticleType,IType>::buildFromPPATagged(int tag,int mask)
{
  console.XDebug() << "TriMesh_PIS_EB::buildFromPPATagged(" << tag << "," << mask << ")\n";
  set<int> id_set;

   
  // for all triangles
  for (
    TriMesh::triangle_iterator tri_iter = this->m_mesh->triangles_begin();
    tri_iter != this->m_mesh->triangles_end();
    tri_iter++
  ){
    // get particles near triangle
    typename ParallelParticleArray<ParticleType>::ParticleListHandle plh=
      ((ParallelParticleArray<ParticleType>*)this->m_ppa)->getParticlesNearTriangle(*tri_iter); 
    console.XDebug() << "triangle " << tri_iter->getID() << " nr. of particles : " << plh->size() << "\n"; 
    // for all particles found
    for(typename ParallelParticleArray<ParticleType>::ParticleListIterator p_iter=plh->begin();
        p_iter!=plh->end();
        p_iter++){
      // if valid interaction
      console.XDebug() << "interaction : " << tri_iter->getID() << " " << (*p_iter)->getID() << "\n";
      if(id_set.find((*p_iter)->getID())==id_set.end()){
        pair<bool,double> dist=tri_iter->dist((*p_iter)->getPos());
        console.XDebug() << "is valid: " << dist.first << " dist : " << dist.second << "\n";
        if(dist.first){
	  int ptag=(*p_iter)->getTag();
	  // if tag is correct, add interaction
          if((ptag & mask)==(tag & mask)){
	    console.XDebug() << "Insert !!\n";
            bool in_flag = this->m_ppa->isInInner((*p_iter)->getPos());
            m_triangle_interactions.push_back(typename IType::TriIntType((*p_iter),&(*tri_iter),m_param,in_flag));
            m_tri_int_set.insert(make_pair(tri_iter->getID(),(*p_iter)->getID()));
            id_set.insert((*p_iter)->getID());
          }
        }
      }
    }
  }
  console.XDebug() << "end  TriMesh_PIS_EB::buildFromPPATagged()";
}

/*!
 */
template<class ParticleType,class IType>
void TriMesh_PIS_EB<ParticleType,IType>::buildFromPPAByGap(double gmax)
{
  console.XDebug() << "TriMesh_PIS_EB::buildFromPPAByGap(" << gmax << ")\n";
  set<int> id_set;

  // for all triangles
  for (
    TriMesh::triangle_iterator tri_iter = this->m_mesh->triangles_begin();
    tri_iter != this->m_mesh->triangles_end();
    tri_iter++
  ){
    // get particles near triangle
    typename ParallelParticleArray<ParticleType>::ParticleListHandle plh=
      ((ParallelParticleArray<ParticleType>*)this->m_ppa)->getParticlesNearTriangle(*tri_iter); 
    console.XDebug() << "triangle " << tri_iter->getID() << " nr. of particles : " << plh->size() << "\n"; 
    // for all particles found
    for(typename ParallelParticleArray<ParticleType>::ParticleListIterator p_iter=plh->begin();
        p_iter!=plh->end();
        p_iter++){
      // if valid interaction
      console.XDebug() << "interaction : " << tri_iter->getID() << " " << (*p_iter)->getID() << "\n";
      if(id_set.find((*p_iter)->getID())==id_set.end()){
        pair<bool,double> dist=tri_iter->dist((*p_iter)->getPos());
        console.XDebug() << "is valid: " << dist.first << " dist : " << dist.second << "\n";
        if(dist.first){
          // check gap
          double gap=fabs(dist.second-(*p_iter)->getRad());
          console.XDebug() << "radius: " << (*p_iter)->getRad() << " gap : " << gap << "\n";
          // if small enough, add
          if(gap<gmax){
            console.XDebug() << "Insert !!\n";
            bool in_flag = this->m_ppa->isInInner((*p_iter)->getPos());
            m_triangle_interactions.push_back(typename IType::TriIntType((*p_iter),&(*tri_iter),m_param,in_flag));
            m_tri_int_set.insert(make_pair(tri_iter->getID(),(*p_iter)->getID()));
            id_set.insert((*p_iter)->getID());
          }
        }
      }
    }
  }
  console.XDebug() << "end  TriMesh_PIS_EB::buildFromPPAByGap()";
}

/*! 
   save snapshot (i.e. viz/postprocess) data

   \param oStream a reference to the stream where the data is written
*/
template<class ParticleType,class IType>
void TriMesh_PIS_EB<ParticleType,IType>::saveSnapShotData(std::ostream &oStream)
{
  const std::string delim = "\n";
  typedef typename IType::TriIntType::CheckPointable CheckPointable;

  // stage 1: count how many interactions with "inner" particles we have
  int icount=0;
  for(typename list<typename IType::TriIntType>::iterator it=m_triangle_interactions.begin();
      it!=m_triangle_interactions.end();
      it++){
    if(it->isInner()) icount++;
  }

  // stage 2: write data
  oStream << IType::getType() << delim;
  oStream << icount << delim;
  for(typename list<typename IType::TriIntType>::iterator it=m_triangle_interactions.begin();
      it!=m_triangle_interactions.end();
      it++){
    if(it->isInner()) CheckPointable(*it).saveCheckPointData(oStream);
    oStream << delim;
  }
}
