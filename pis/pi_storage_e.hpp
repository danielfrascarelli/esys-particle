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

// STL includes
#include <algorithm>  // for sort, copy
#include <iterator> // for back_inserter

using std::sort;
using std::copy;
using std::back_inserter;

template<typename P,typename InteractionType>
const int ParallelInteractionStorage_E<P,InteractionType>::m_exchg_tag=43;

/*!
  Construct parallel interaction storage

  \param PPA a pointer to the particle array
  \param param the interaction parameters
*/
template<typename P,typename I>
ParallelInteractionStorage_E<P,I>::ParallelInteractionStorage_E(
  AParallelParticleArray* PPA,
  const typename I::ParameterType& param
)
  : TParallelInteractionStorage<I>(PPA),
    m_comm(PPA->getComm()),
    m_param(param)
{
  m_unbreakable=false;
}

/*!

*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_E<P,InteractionType>::exchange()
{
  for(int i=0;i<3;i++){
    if(m_comm.get_dim(i)>1){
      // -- up --
      exchange_boundary(i,1);
      // -- down --
      exchange_boundary(i,-1);
    }
  } 
}

/*!
  helper function to do the actual shifting of values in exchange()

  \param dim dimension, 0->x, 1->y, 2->z
  \param dir direction, 1->up, -1->down
*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_E<P,InteractionType>::exchange_boundary(int dim,int dir)
{ 
  console.XDebug() << "PIS_E::exchange_boundary(" << dim << "," << dir << ") at node " << m_comm.rank() << "\n";
  set<int> bdry_ids;
  vector<InteractionType> recv_buffer;
  vector<InteractionType> send_buffer;

  // get boundary
  bdry_ids = this->m_ppa->getBoundarySlabIds(dim,dir);
  // for all interactions
  for(
      typename list<InteractionType>::iterator iter = this->m_interactions.begin();
      iter != this->m_interactions.end();
      iter++
  ){
    vector<int> pids=iter->getAllID(); // get particle IDs
    bool flag=false;
    // check if any id is in boundary slab
    vector<int>::iterator it2=pids.begin();
    while(it2!=pids.end() && !flag){
      flag=(bdry_ids.find(*it2)!=bdry_ids.end());
      it2++;
    }
    if(flag){
      send_buffer.push_back(*iter);
    }
  }
  // shift 
  m_comm.shift_cont_packed(send_buffer,recv_buffer,dim,dir,m_exchg_tag);
  // try to insert the received interactions
  for(typename vector<InteractionType>::iterator iter=recv_buffer.begin();
      iter!=recv_buffer.end();
      iter++){
    tryInsert(*iter);
  }
  // clean buffers
  send_buffer.clear();
  recv_buffer.clear();
  console.XDebug() << "end PIS_E::exchange_boundary(" << dim << "," << dir << ") at node " << m_comm.rank() << "\n";
}

/*!
  rebuild all interactions, i.e. set particle pointers according to 
  particle indices
*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_E<P,InteractionType>::rebuild()
{
  console.XDebug() << "PIS_E::rebuild at node " << m_comm.rank() << "\n";
  console.XDebug() << "size pre rebuild: " << this->m_interactions.size() << "\n";

  // -- DEBUG ---
  for(typename list<InteractionType>::iterator iter = this->m_interactions.begin();
      iter!=this->m_interactions.end();
      iter++){
    vector<int> pids=iter->getAllID();
    console.XDebug() << pids[0] << " - " << pids[1] << "\n";
  }
  // --- END DEBUG ---
  vector<P*> pptr;
  ParallelParticleArray<P>* t_ppa=(ParallelParticleArray<P>*)(this->m_ppa);
  typename list<InteractionType>::iterator iter = this->m_interactions.begin();
  while(iter != this->m_interactions.end()){
    vector<int> pids=iter->getAllID();
    vector<int>::const_iterator it2=pids.begin();
    bool flag=true;
    // check if the particles with the stored IDs are here
    while(it2!=pids.end() && flag){ 
      P* ptr=t_ppa->getParticlePtrByIndex(*it2);
      if(ptr!=NULL){
	pptr.push_back(ptr);
      } else {
	flag=false;
      }
      it2++;
    }
    if(flag){ // if all particle IDs are valid -> set particle pointers 
      iter->setPP(pptr);
      iter->checkIDs();
      iter++;
    } else { // if not -> erase interactions
      const typename list<InteractionType>::iterator er_iter=iter;
      iter++;
      this->m_interactions.erase(er_iter);
      m_set.erase(make_pair(pids[0],pids[1]));
    }
    pptr.clear();
  }
  console.XDebug() << "size post rebuild: " << this->m_interactions.size() << "\n";
  //  cout << "end PIS_E::rebuild at node " << m_comm.rank() << endl;
}


/*!
  Try to insert interaction into the PIS. If insertion is impossible because the interaction 
  is already in, or one of the particles is not in the associated PPA nothing happens.

  \param In the interaction
*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_E<P,InteractionType>::tryInsert(const InteractionType& In)
{
  bool flag=true;
  
  ParallelParticleArray<P>* t_ppa=(ParallelParticleArray<P>*)(this->m_ppa);
  // check if interaction is already in 
  vector<int> pids=In.getAllID();
  flag=!isIn(pids);
  // try to get particle pointers from ppa 
  vector<int>::const_iterator iter=pids.begin();
  while(iter!=pids.end() && flag){
     P* ptr=t_ppa->getParticlePtrByIndex(*iter);
     if(ptr!=NULL){
       //pptr.push_back(ptr);
    } else {
      flag=false;
    }
    iter++;
  }
 
  if(flag){
    this->m_interactions.push_back(In);
    m_set.insert(make_pair(pids[0],pids[1]));
  } 
}

/*!
  \warning evil hack, only checks 1st & 2nd id -> change from pair<int,int> to vector<int>
*/
template<typename P,typename InteractionType>
bool ParallelInteractionStorage_E<P,InteractionType>::isIn(const vector<int>& pids)
{
  bool res;

  if(pids[0] > pids [1]){
    console.Debug()<< "flipped PIDS : " << pids[0] << "," << pids[1] << "\n";
  }
  res=m_set.find(make_pair(pids[0],pids[1]))!=m_set.end();

  return res;
}


/*!
  Insert interactions newly created from particle Ids and parameters. If insertion is impossible 
  because the interaction is already in, or one of the particles is not in the associated PPA 
  nothing happens.

  \param cpids the particle Ids
*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_E<P,InteractionType>::tryInsert(const vector<int>& cpids)
{
  vector<P*> pptr;
  bool flag=true;
  
  // can't sort pids directly because of const -> need copy
  vector<int> pids;
  copy(cpids.begin(),cpids.end(),back_inserter(pids));
  sort(pids.begin(),pids.end());
  
  ParallelParticleArray<P>* t_ppa=(ParallelParticleArray<P>*)(this->m_ppa);
  // check if interaction is already in 
  flag=!isIn(pids);
  // try to get particle pointers from ppa 
  vector<int>::const_iterator iter=pids.begin();
  while(iter!=pids.end() && flag){
     P* ptr=t_ppa->getParticlePtrByIndex(*iter);
     if(ptr!=NULL){
      pptr.push_back(ptr);
    } else {
      flag=false;
    }
    iter++;
  }
 
  if(flag){
    // initialize interaction from particle pointers and interaction parameters
    InteractionType new_interaction(pptr[0],pptr[1],m_param);
    vector<int> allid=new_interaction.getAllID();
    console.XDebug() << allid[0] << " , " << allid[1] << "\n"; 
    // insert interaction
    this->m_interactions.push_back(new_interaction);
    this->m_set.insert(make_pair(pids[0],pids[1]));
  } 
}

/*!
  calculate all forces
*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_E<P,InteractionType>::calcForces()
{
  console.Debug()
    << "calculating "
    << this->m_interactions.size()
    << " interaction forces\n" ;
  
  for(
    typename list<InteractionType>::iterator it = this->m_interactions.begin();
    it != this->m_interactions.end();
    it++
  ){
    it->calcForces();
  }
}

/*!
    set the interactions "unbreakable" -> turns update into a NO-OP
 
      \param b true -> unbreakable, false -> breakable
*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_E<P,InteractionType>::setUnbreakable(bool b)
{
  m_unbreakable=b;
}
