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

template<class ParticleType,class IType>
const int Mesh2D_PIS_EB<ParticleType,IType>::m_edge_exchg_tag=45;
template<class ParticleType,class IType>
const int Mesh2D_PIS_EB<ParticleType,IType>::m_corner_exchg_tag=46;

/*!
  constructor

  \param mesh_p a pointer to the triangle mesh
  \param ppa_p a pointer to the particle array
  \param param the interaction parameters
*/
template<class ParticleType,class IType>
Mesh2D_PIS_EB<ParticleType,IType>::Mesh2D_PIS_EB(Mesh2D* mesh_p,ParallelParticleArray<ParticleType>* ppa_p,typename IType::ParameterType param)
  :Mesh2D_PIS<ParticleType>(mesh_p,ppa_p),m_comm(ppa_p->getComm())
{
  console.XDebug() << "Mesh2D_PIS_EB constructor\n";
  m_param=param;
  this->m_update_timestamp=0;
}
   
/*!
  Check if an interaction is in this PIS. The first 2 values in the vector
  are expected to be the edge/corner (v[0]) and particle (v[1]) ids, the 3rd 
  an indicator if tri edge (v[2]==1)or corner (v[2]==2) interaction. 
  If there is no 3rd value or it is not 1 (edge), "false" is returned.

  \param v vector of particle ids
  \warning log(N)
*/
template<class ParticleType,class IType>
bool Mesh2D_PIS_EB<ParticleType,IType>::isIn(const std::vector<int>& v)
{
  bool res=false;
  
  if(v.size()<3){
    res=false;
  } else {
    switch (v[2]){
    case 1: res=m_edge_int_set.find(make_pair(v[0],v[1]))!=m_edge_int_set.end(); break;
    case 2: res=m_edge_int_set.find(make_pair(v[0],v[1]))!=m_corner_int_set.end(); break;
    default: console.Error() << "wrong value in argument of Mesh2D_PIS::isIn !!\n"; break;
    }
  }

  return res;
}

/*!
  calculate all the forces
*/
template<class ParticleType,class IType>
void Mesh2D_PIS_EB<ParticleType,IType>::calcForces()
{
  console.XDebug() << "Mesh2D_PIS_EB calculating " << m_edge_interactions.size() << " edge forces and"
		   << m_corner_interactions.size() << " corner forces\n";

  // calculate forces for edge interactions
  for(typename list<typename IType::TriIntType>::iterator ed_iter=m_edge_interactions.begin();
      ed_iter!=m_edge_interactions.end();
      ed_iter++){
    ed_iter->calcForces();
  }
  // calculate forces for corner interactions
  for(typename list<typename IType::CornerIntType>::iterator c_iter=m_corner_interactions.begin();
      c_iter!=m_corner_interactions.end();
      c_iter++){
    c_iter->calcForces();
  }
}

/*!
 */
template<class ParticleType,class IType>
bool Mesh2D_PIS_EB<ParticleType,IType>::update()
{
  console.XDebug() << "Mesh2D_PIS_EB::update on node " << m_comm.rank() << "\n"; 
  bool res=false;

  // edge interactions
  typename list<typename IType::TriIntType>::iterator iter=m_edge_interactions.begin();
  while(iter!=m_edge_interactions.end()){
    if(iter->broken()){
      res=true;
      typename list<typename IType::TriIntType>::iterator er_iter=iter;
      iter++;
      // remove ids from map
      m_edge_int_set.erase(make_pair(er_iter->getTid(),er_iter->getPid()));
      // remove interaction
      m_edge_interactions.erase(er_iter);
    } else {
      iter++;
    }
  }
  // corner interactions
  typename list<typename IType::CornerIntType>::iterator c_iter=m_corner_interactions.begin();
  while(c_iter!=m_corner_interactions.end()){
    if(c_iter->broken()){
      res=true;
      typename list<typename IType::CornerIntType>::iterator cr_iter=c_iter;
      c_iter++;
      // remove ids from map
      m_corner_int_set.erase(make_pair(cr_iter->getCid(),cr_iter->getPid()));
      // remove interaction
      m_corner_interactions.erase(cr_iter);
    } else {
      c_iter++;
    }
  }
  console.XDebug() << "end Mesh2D_PIS_EB::update on node " << m_comm.rank() << "\n"; 
  return res;
}


/*!
 */
template<class ParticleType,class IType>
void Mesh2D_PIS_EB<ParticleType,IType>::exchange()
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
  helper function to do the actual shifting of values in exchange()

  \param dim dimension, 0->x, 1->y, 2->z
  \param dir direction, 1->up, -1->down
*/
template<class ParticleType,class IType>
void Mesh2D_PIS_EB<ParticleType,IType>::exchange_boundary(int dim,int dir)
{
  console.XDebug() << "Mesh2D_PIS_EB::exchange_boundary(" << dim << "," << dir << ") at node " << m_comm.rank() << "\n";
  
  std::set<int> bdry_ids;
  std::vector<typename IType::TriIntType> recv_tri_buffer;
  std::vector<typename IType::TriIntType> send_tri_buffer;
  std::vector<typename IType::CornerIntType> recv_corner_buffer;
  std::vector<typename IType::CornerIntType> send_corner_buffer;

  // --- setup data to send ---
  // get boundary
  bdry_ids = this->m_ppa->getBoundarySlabIds(dim,dir);
  // --- edges ---
  // for all interactions
  for(typename std::list<typename IType::TriIntType>::iterator iter=m_edge_interactions.begin();
      iter!=m_edge_interactions.end();
      iter++){
    int pid=iter->getPid(); 
    if(bdry_ids.find(pid)!=bdry_ids.end()) { // if particle in interaction is in bdry -> put in buffer
      send_tri_buffer.push_back(*iter);
    }  
  }
  // shift 
  m_comm.shift_cont_packed(send_tri_buffer,recv_tri_buffer,dim,dir,m_edge_exchg_tag);
  // insert received data
  for(typename std::vector<typename IType::TriIntType>::iterator iter=recv_tri_buffer.begin();
      iter!=recv_tri_buffer.end();
      iter++){
    tryInsert(*iter);
  }
  // --- corners ---
  // for all interactions
  for(typename std::list<typename IType::CornerIntType>::iterator iter=m_corner_interactions.begin();
       iter!=m_corner_interactions.end();
       iter++){
    int pid=iter->getPid(); 
    if(bdry_ids.find(pid)!=bdry_ids.end()) { // if particle in interaction is in bdry -> put in buffer
      send_corner_buffer.push_back(*iter);
    }  
  }
  // shift 
  m_comm.shift_cont_packed(send_corner_buffer,recv_corner_buffer,dim,dir,m_corner_exchg_tag);
  // insert received data
  for(typename std::vector<typename IType::CornerIntType>::iterator iter=recv_corner_buffer.begin();
      iter!=recv_corner_buffer.end();
      iter++){
    tryInsert(*iter);
  }
  
  console.XDebug() << "end Mesh2D_PIS_EB::exchange_boundary\n";
}

template<class ParticleType,class IType>
void Mesh2D_PIS_EB<ParticleType,IType>::setTimeStepSize(double dt)
{
}

/*!
  Rebuild interactions after moving particles or interactions between processes. Set 
  particle pointers accordig to particle IDs and remove interactionw which include
  unavailable particles.
*/
template<class ParticleType,class IType>
void Mesh2D_PIS_EB<ParticleType,IType>::rebuild()
{
  console.XDebug() << "Mesh2D_PIS_EB::rebuild on node " << m_comm.rank() << "\n"; 
  ParallelParticleArray<ParticleType>* t_ppa =
    (ParallelParticleArray<ParticleType>*)(this->m_ppa); // should be a dynamic_cast
  // for all edge interactions
  typename std::list<typename IType::TriIntType>::iterator ti_iter=m_edge_interactions.begin();
  while(ti_iter!=m_edge_interactions.end()){
    int pid=ti_iter->getPid();
    ParticleType *part_p=t_ppa->getParticlePtrByIndex(pid);
    if(part_p!=NULL) { // particle is available in m_ppa -> setup pointer in interaction
      ti_iter->setPP(part_p);
      Edge2D *ed_p = this->m_mesh->getEdgeById(ti_iter->getTid());
      ti_iter->setTP(ed_p);
      ti_iter++;
    } else { // particle is not available in m_ppa -> remove interaction
      const typename list<typename IType::TriIntType>::iterator er_iter=ti_iter;
      ti_iter++;
      m_edge_int_set.erase(make_pair(er_iter->getTid(),pid));
      m_edge_interactions.erase(er_iter); 
    }
  }
  // and now for the corners
  typename list<typename IType::CornerIntType>::iterator ci_iter=m_corner_interactions.begin();
  while(ci_iter!=m_corner_interactions.end()){
    int pid=ci_iter->getPid();
    ParticleType *part_p=t_ppa->getParticlePtrByIndex(pid);
    if(part_p!=NULL) { // particle is available in m_ppa -> setup pointer in interaction
      ci_iter->setPP(part_p);
      Corner2D *co_p = this->m_mesh->getCornerById(ci_iter->getCid());
      ci_iter->setCP(co_p);
      ci_iter++;
    } else { // particle is not available in m_ppa -> remove interaction
      const typename list<typename IType::CornerIntType>::iterator cr_iter=ci_iter;
      ci_iter++;
      m_corner_int_set.erase(make_pair(cr_iter->getCid(),pid));
      m_corner_interactions.erase(cr_iter); 
    }
  }

  console.XDebug() << "end Mesh2D_PIS_EB::rebuild on node " << m_comm.rank() << "\n"; 
}
   
/*!
 */
template<class ParticleType,class IType>
void Mesh2D_PIS_EB<ParticleType,IType>::tryInsert(const typename IType::TriIntType& In)
{
  console.XDebug() << "Mesh2D_PIS_EB::tryInsert(const typename IType::TriIntType& In)\n";
  ParallelParticleArray<ParticleType>* t_ppa =
    (ParallelParticleArray<ParticleType>*)(this->m_ppa);
  // check if interaction is already in 
  bool is_in=(m_edge_int_set.find(make_pair(In.getTid(),In.getPid()))!=m_edge_int_set.end());
  ParticleType *part_p=t_ppa->getParticlePtrByIndex(In.getPid());
  if((!is_in) && (part_p!=NULL)){
    m_edge_interactions.push_back(In);
    m_edge_int_set.insert(make_pair(In.getTid(),In.getPid()));
  }
}

/*!
 */
template<class ParticleType,class IType>
void Mesh2D_PIS_EB<ParticleType,IType>::tryInsert(const typename IType::CornerIntType& In)
{
  console.XDebug() << "Mesh2D_PIS_EB::tryInsert(const typename IType::CornerIntType& In)\n";
  ParallelParticleArray<ParticleType>* t_ppa =
    (ParallelParticleArray<ParticleType>*)(this->m_ppa);
  // check if interaction is already in 
  bool is_in=(m_corner_int_set.find(make_pair(In.getCid(),In.getPid()))!=m_corner_int_set.end());
  ParticleType *part_p=t_ppa->getParticlePtrByIndex(In.getPid());
  if((!is_in) && (part_p!=NULL)){
    m_corner_interactions.push_back(In);
    m_corner_int_set.insert(make_pair(In.getCid(),In.getPid()));
  }
}

/*!
  Insert interactions newly created from particle Ids and parameters. If insertion is impossible 
  because the interaction is already in, or one of the particles is not in the associated PPA 
  nothing happens.
  Check if an interaction is in this PIS. The first 2 values in the vector
  are expected to be the tri/edge/corner (pids[0]) and particle (pids[1]) ids, the 3rd 
  an indicator if edge (pids[2]==1)or corner (pids[2]==2) interaction. 
  If there is no 3rd value or it is not in [1,2], nothing happens.

  \param pids the particle Ids
*/
template<class ParticleType,class IType>
void Mesh2D_PIS_EB<ParticleType,IType>::tryInsert(const vector<int>& pids)
{
  console.XDebug() << "Mesh2D_PIS_EB::(const vector<int>& pids)\n";
  ParallelParticleArray<ParticleType>* t_ppa =
    (ParallelParticleArray<ParticleType>*)(this->m_ppa); // should be dynamic_cast

  if(pids.size()<3){
    bool is_in=isIn(pids); // interaction already in
    ParticleType *part_p=t_ppa->getParticlePtrByIndex(pids[1]);
    if((!is_in) && (part_p!=NULL)){ 
      switch (pids[2]){
      case 1: { // edge
        Edge2D *edge_p = this->m_mesh->getEdgeById(pids[0]);
        if(edge_p!=NULL){
          m_edge_interactions.push_back(
            typename IType::TriIntType(
              part_p,
              edge_p,
              m_param,
              this->m_ppa->isInInner(part_p->getPos())
            )
          );
          m_edge_int_set.insert(make_pair(pids[0],pids[1]));
        } else {
          console.Error() << "ERROR: Wrong edge id " << pids[0] << " in Mesh2D_PIS_EB::tryInsert\n";
        }
      } break;
      case 2: {
        Corner2D *corner_p = this->m_mesh->getCornerById(pids[0]);
        if(corner_p!=NULL){
          m_corner_interactions.push_back(
            typename IType::CornerIntType(
              part_p,
              corner_p,
              m_param,
              this->m_ppa->isInInner(part_p->getPos())
            )
          );
          m_corner_int_set.insert(make_pair(pids[0],pids[1]));
        } else {
          console.Error() << "ERROR: Wrong corner id " << pids[0] << " in Mesh2D_PIS_EB::tryInsert\n";
        }
      } break;
      default : console.Error() << "Error: pids[2]= " << pids[2] << "\n"; // Can't happen !! 
      }
    }
  }
}


/*!
 */
template<class ParticleType,class IType>
void Mesh2D_PIS_EB<ParticleType,IType>::buildFromPPATagged(int tag ,int mask)
{
  console.XDebug() << "Mesh2D_PIS_EB::buildFromPPATagged(" << tag << "," << mask << ")\n";
  set<int> id_set;  

  // for all edges
  for(
    Mesh2D::edge_iterator ed_iter = this->m_mesh->edges_begin();
    ed_iter != this->m_mesh->edges_end();
    ed_iter++){
    // get particles near edge
    typename ParallelParticleArray<ParticleType>::ParticleListHandle plh=
      ((ParallelParticleArray<ParticleType>*)this->m_ppa)->getParticlesNearEdge(&(*ed_iter)); 
    console.Debug() << "edge " << ed_iter->getID() << " nr. of particles : " << plh->size() << "\n"; 
    // for all particles found
    for(
      typename ParallelParticleArray<ParticleType>::ParticleListIterator p_iter=plh->begin();
      p_iter!=plh->end();
      p_iter++){
      // if valid interaction
      console.Debug() << "interaction : " << ed_iter->getID() << " " << (*p_iter)->getID() << "\n";
      if(id_set.find((*p_iter)->getID())==id_set.end()){
        pair<bool,double> dist=ed_iter->dist((*p_iter)->getPos());
        console.Debug() << "is valid: " << dist.first << " dist : " << dist.second << "\n";
        if(dist.first){
	  int ptag=(*p_iter)->getTag();
	  // if tag is correct, add interaction
          if((ptag & mask)==(tag & mask)){
            console.Debug() << "Inserting " << (*p_iter)->getID() << " !!\n";
            bool in_flag = this->m_ppa->isInInner((*p_iter)->getPos());
            m_edge_interactions.push_back(typename IType::TriIntType((*p_iter),&(*ed_iter),m_param,in_flag));
            m_edge_int_set.insert(make_pair(ed_iter->getID(),(*p_iter)->getID()));
            id_set.insert((*p_iter)->getID());
          }
        }
      }
    }
    
  }
  console.XDebug() << "end  Mesh2D_PIS_EB::buildFromPPATagged()";
}

/*!
  build interactions according to given maximum gap between particle and 2d edge

  \param gmax the maximum gap
*/
template<class ParticleType,class IType>
void Mesh2D_PIS_EB<ParticleType,IType>::buildFromPPAByGap(double gmax)
{
  console.XDebug() << "Mesh2D_PIS_EB::buildFromPPAByGap(" << gmax << ")\n";
  set<int> id_set;

  // for all edges
  for(
    Mesh2D::edge_iterator ed_iter = this->m_mesh->edges_begin();
    ed_iter != this->m_mesh->edges_end();
    ed_iter++
  ){
    // get particles near edge
    typename ParallelParticleArray<ParticleType>::ParticleListHandle plh=
      ((ParallelParticleArray<ParticleType>*)this->m_ppa)->getParticlesNearEdge(&(*ed_iter)); 
    console.Debug() << "edge " << ed_iter->getID() << " nr. of particles : " << plh->size() << "\n"; 
    // for all particles found
    for(
      typename ParallelParticleArray<ParticleType>::ParticleListIterator p_iter=plh->begin();
      p_iter!=plh->end();
      p_iter++){
      // if valid interaction
      console.Debug() << "interaction : " << ed_iter->getID() << " " << (*p_iter)->getID() << "\n";
      if(id_set.find((*p_iter)->getID())==id_set.end()){
        pair<bool,double> dist=ed_iter->dist((*p_iter)->getPos());
        console.Debug() << "is valid: " << dist.first << " dist : " << dist.second << "\n";
        if(dist.first){
          // check ga
          double gap=fabs(dist.second-(*p_iter)->getRad());
          console.Debug() << "radius: " << (*p_iter)->getRad() << " gap : " << gap << "\n";
          // if small enough, add
          if(gap<gmax){
            console.Debug() << "Inserting " << (*p_iter)->getID() << " !!\n";
            bool in_flag = this->m_ppa->isInInner((*p_iter)->getPos());
            m_edge_interactions.push_back(typename IType::TriIntType((*p_iter),&(*ed_iter),m_param,in_flag));
            m_edge_int_set.insert(make_pair(ed_iter->getID(),(*p_iter)->getID()));
            id_set.insert((*p_iter)->getID());
          }
        }
      }
    }
  }
}

template <class ParticleType,class IType> 
typename Mesh2D_PIS_EB<ParticleType,IType>::InteractionIterator  Mesh2D_PIS_EB<ParticleType,IType>::getInnerInteractionIterator()
{
  return 
    InteractionIterator(
      m_edge_interactions.begin(),
      m_edge_interactions.end(),
      this->m_ppa
  );
}

template <class ParticleType,class IType> 
void Mesh2D_PIS_EB<ParticleType,IType>::saveSnapShotData(std::ostream& ost)
{
  const std::string delim = "\n";
  typedef typename IType::TriIntType::CheckPointable CheckPointable;

  InteractionIterator it = getInnerInteractionIterator();
  ost << IType::getType() << delim;
  ost << it.getNumRemaining();
  while (it.hasNext()){
    ost << delim;
    CheckPointable(it.next()).saveSnapShotData(ost);
  }
}

template <class ParticleType,class IType>
void Mesh2D_PIS_EB<ParticleType,IType>::saveCheckPointData(std::ostream& ost)
{
  const std::string delim = "\n";
  typedef typename IType::TriIntType::CheckPointable CheckPointable;

  InteractionIterator it = getInnerInteractionIterator();
  ost << IType::getType() << delim;
  ost << it.getNumRemaining();
  while (it.hasNext()){
    ost << delim;
    CheckPointable(it.next()).saveCheckPointData(ost);
  }
}

template <class ParticleType,class IType> 
void Mesh2D_PIS_EB<ParticleType,IType>::loadCheckPointData(std::istream& ost)
{
  console.Critical() << "Mesh2D_PIS_EB<ParticleType,IType>::loadCheckPointData NOT IMPLEMENTED\n";
}

// --- Iterator members ---

template <class ParticleType,class IType> 
Mesh2D_PIS_EB<ParticleType,IType>::InteractionIterator::InteractionIterator(Iterator begin, Iterator end, AParallelParticleArray* ppa)
  : m_numRemaining(0),
    m_it(end),
    m_end(end),
    m_ppa(ppa)
{
  m_numRemaining = 0;
  for (Iterator it = begin; it != end; it++) {
    if  (isInner(it)) {
      m_numRemaining++;
    }
  }
  m_it  = begin;
  m_end = end;
}


template <class ParticleType,class IType> 
bool Mesh2D_PIS_EB<ParticleType,IType>::InteractionIterator::hasNext()
{
  return (m_numRemaining > 0);
}

template <class ParticleType,class IType> 
typename IType::TriIntType& Mesh2D_PIS_EB<ParticleType,IType>::InteractionIterator::next()
{
  while (!isInner(m_it)) {
    m_it++;
  }
  Interaction &i = *m_it;
  m_it++;
  m_numRemaining--;
  return i;
}

template <class ParticleType,class IType> 
int Mesh2D_PIS_EB<ParticleType,IType>::InteractionIterator::getNumRemaining()
{
  return m_numRemaining;
}

template <class ParticleType,class IType> 
bool Mesh2D_PIS_EB<ParticleType,IType>::InteractionIterator::isInner(const Iterator& it)
{
  return m_ppa->isInInner(it->getPos());
}
