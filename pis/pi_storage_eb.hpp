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

template<typename P,typename I>
ParallelInteractionStorage_EB<P,I>::ParallelInteractionStorage_EB(AParallelParticleArray* ppa,const typename I::ParameterType& param):ParallelInteractionStorage_E<P,I>(ppa,param)
{
  m_unbreakable=false;
}

/*!
  Update interactions. Check for broken interactions and remove them.
*/
template<typename P,typename InteractionType>
bool ParallelInteractionStorage_EB<P,InteractionType>::update()
{
  bool res=false;

  if(!m_unbreakable){
    console.XDebug() << "PIS_E::updating\n"; 
    typename list<InteractionType>::iterator iter = this->m_interactions.begin();
    while (iter != this->m_interactions.end()){
      if(iter->broken()){
	res=true;
	typename list<InteractionType>::iterator er_iter=iter;
	// get IDs to remove from set
	vector<int> pids=iter->getAllID();
	this->m_set.erase(make_pair(pids[0],pids[1]));
	iter++;
	// remove interaction
	this->m_interactions.erase(er_iter);
      } else {
	iter++;
      }
    }
  } else {
    console.XDebug() << "PIS_E::not updating\n"; 
  }

  return res;
}

/*! 
   save snapshot (i.e. viz/postprocess) data
*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_EB<P,InteractionType>::saveSnapShotData(std::ostream &oStream)
{
  const std::string delim = "\n";
  typedef typename InteractionType::CheckPointable CheckPointable;

  typename ParallelInteractionStorage_E<P,InteractionType>::InteractionIterator it =
    this->getInnerInteractionIterator();
  oStream << InteractionType::getType() << delim;
  oStream << it.getNumRemaining();
  if (it.hasNext()) {
    oStream << delim;
    CheckPointable(it.next()).saveCheckPointData(oStream);
    while (it.hasNext())
    {
      oStream << delim;
      CheckPointable(it.next()).saveCheckPointData(oStream);
    }
  }
}

/*! 
   save checkpoint (i.e. restart) data
*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_EB<P,InteractionType>::saveCheckPointData(std::ostream &oStream)
{
  const std::string delim = "\n";
  //  typedef typename InteractionType::CheckPointable CheckPointable;

  typename ParallelInteractionStorage_E<P,InteractionType>::InteractionIterator it =
    this->getInnerInteractionIterator();
  oStream << InteractionType::getType() << delim;
  oStream << it.getNumRemaining();
  if (it.hasNext()) {
    oStream << delim;
    it.next().saveRestartData(oStream);
    while (it.hasNext())
    {
      oStream << delim;
      it.next().saveRestartData(oStream);
    }
  }
}

/*!
  Read interaction data from input stream pointing to a restartable checkpoint file.
  The stream needs to be already positioned at the right place.
 
  \param iStream the input stream
  \warning return type may change to bool at some stage
*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_EB<P,InteractionType>::loadCheckPointData(std::istream &iStream)
{
  // read interaction type from stream
  std::string cp_interaction_type;
  iStream >> cp_interaction_type;
  // compare interaction type in stream with type of this IG
  // in not equal, signal error 
  if(cp_interaction_type!=InteractionType::getType()){
    std::cerr << "interaction types differ between checkpoint " 
	      << cp_interaction_type << " and scipt "
	      << InteractionType::getType() << std::endl;
  } else { // correct type -> read data
    // read nr. of bonds in IG
    int nconn;
    iStream >> nconn;
    std::cerr << "reading " << nconn << "  " << InteractionType::getType() << " interactions " << std::endl;

    ParallelParticleArray<P>* t_ppa=(ParallelParticleArray<P>*)(this->m_ppa);
    
    // -- read bonds
    for(int i=0;i<nconn;i++){
      InteractionType new_bond;
      // read a bond
      new_bond.loadRestartData(iStream);
      // set particle pointers
      vector<int> pids=new_bond.getAllID();
      P* ptr1=t_ppa->getParticlePtrByIndex(pids[0]);
      P* ptr2=t_ppa->getParticlePtrByIndex(pids[1]);
      if((ptr1!=NULL) && (ptr2!=NULL)){
	new_bond.setPP(ptr1,ptr2);
      } else {
	std::cerr << "trying to insert bond: particles with Id " << pids[0] << " , " << pids[1] << "not present!" << std::endl;
      }
      // insert it into interaction storage
      this->tryInsert(new_bond);
    }
  }
}

template<typename P,typename InteractionType>
void ParallelInteractionStorage_EB<P,InteractionType>::calcHeatTrans()
{
  console.Debug()
    << "calculating " << this->m_interactions.size()
    << " heat interaction transfers\n" ;

  for(
    typename list<InteractionType>::iterator it = this->m_interactions.begin();
    it != this->m_interactions.end();
    it++
  ){
    it->calcHeatTrans();
  }
}

/*!
  set the interactions "unbreakable" -> turns update into a NO-OP

  \param b true -> unbreakable, false -> breakable
*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_EB<P,InteractionType>::setUnbreakable(bool b)
{
  m_unbreakable=b;
}
