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
ParallelInteractionStorage_ED<P,I>::ParallelInteractionStorage_ED(AParallelParticleArray* ppa,const typename I::ParameterType& param):ParallelInteractionStorage_E<P,I>(ppa,param)
{
  m_update_timestamp=0;
}

template<typename P,typename InteractionType>
void  ParallelInteractionStorage_ED<P,InteractionType>::addExIG(AParallelInteractionStorage* eg)
{
  console.Debug() << "setExIG " << eg << "\n";
  m_exIG.push_back(eg); 
  // clean out excluded interactions
  typename list<InteractionType>::iterator iter = this->m_interactions.begin();
  while(iter != this->m_interactions.end()){
    // check if in exIG
    vector<int> rm_pids=iter->getAllID();
    if(eg->isIn(rm_pids)){
      console.XDebug() << "removing excluded: " << rm_pids[0] << " - " << rm_pids[1] << "\n";
      typename list<InteractionType>::iterator er_iter=iter;
      // get particle ids and remove pair from set
      this->m_set.erase(make_pair(rm_pids[0],rm_pids[1]));
      iter++;
      this->m_interactions.erase(er_iter);
    } else {
      iter++;
    }
  }
}


/*!
  check if a potential interaction is in one of the excluding IGs

  \param tv the potential interaction
*/
template<typename T,typename InteractionType>
bool ParallelInteractionStorage_ED<T,InteractionType>:: isExcluded(const vector<int> tv)
{ 
  bool in_exig=false;
  if(m_exIG.size()>0){ // if there is an ExIG
    vector<AParallelInteractionStorage*>::iterator exiter=m_exIG.begin();
    while(exiter!=m_exIG.end() && (!in_exig)){ // not yet excluded and more excluders to check 
      in_exig=(*exiter)->isIn(tv);
      exiter++;
    }
  }
  return in_exig;
}

template<typename T,typename InteractionType>
void ParallelInteractionStorage_ED<T,InteractionType>::setTimeStepSize(
  double dt
)
{
  console.Debug()
    << "setting time step size for "
    << this->m_interactions.size() << " interaction forces\n" ;
  this->m_param.setTimeStepSize(dt);
  for (
    typename std::list<InteractionType>::iterator it = this->m_interactions.begin();
    it != this->m_interactions.end();
    it++
  ){
    it->setTimeStepSize(dt);
  }
}

/*!
  Update interactions. Do full dynamic search.
*/
template<typename T,typename InteractionType>
bool ParallelInteractionStorage_ED<T,InteractionType>::update()
{
  //std::cout << "ParallelInteractionStorage_ED::update at node " << m_comm.rank() << std::endl << std::flush;
  int count_l=0;
  bool res=true;

  if (this->m_update_timestamp != this->m_ppa->getTimeStamp()){// m_ppa rebuild since last update 
    console.XDebug() << "node " << this->m_comm.rank() << " ppa has been rebuilt\n";
    // clean out old interactions if not flagged as persistent
    typename list<InteractionType>::iterator iter = this->m_interactions.begin();
    while(iter != this->m_interactions.end()){
      if(iter->isPersistent()){
        iter++;
        //console.XDebug() << "node " << m_comm.rank() << "persistent interaction\n";
      }else{
        typename list<InteractionType>::iterator er_iter=iter;
        // get particle ids and remove pair from set
        vector<int> rm_pids=iter->getAllID();
        this->m_set.erase(make_pair(rm_pids[0],rm_pids[1]));
        iter++;
        this->m_interactions.erase(er_iter);
      }
    }
    // get list  of pairs from m_ppa
    typename ParallelParticleArray<T>::PairListHandle plh =
      ((ParallelParticleArray<T>*)this->m_ppa)->getFullPairList();
    // generate interactions from pairs
    for(typename ParallelParticleArray<T>::PairListIterator iter=plh->begin();
	iter!=plh->end();
	iter++){
      // check vs. ExIG
      vector<int> tv;
      // ids in pair
      int id1=iter->first->getID();
      int id2=iter->second->getID();
      tv.push_back(id1);
      tv.push_back(id2);
      if((!isExcluded(tv))&&(!this->isIn(tv))){  // if not already in or in ExIG
	this->m_interactions.push_back(
	  InteractionType(iter->first,iter->second,this->m_param)
	  );
	this->m_set.insert(make_pair(id1,id2));
	count_l++;
      }
    }
  } else { // m_ppa not rebuild since last update -> just get additional interactions
    console.XDebug() << "node " << this->m_comm.rank() << " ppa not rebuilt\n";
    // get list  of pairs from m_ppa
    typename ParallelParticleArray<T>::PairListHandle plh =
      ((ParallelParticleArray<T>*)this->m_ppa)->getNewPairList();
    for (
      typename ParallelParticleArray<T>::PairListIterator iter=plh->begin();
      iter!=plh->end();
      iter++
    ){
      // check vs. ExIG
      vector<int> tv;
      // ids in pair
      int id1=iter->first->getID();
      int id2=iter->second->getID();
      tv.push_back(id1);
      tv.push_back(id2); 
      if((!isExcluded(tv))&&(!this->isIn(tv))){
	this->m_interactions.push_back(
            InteractionType(iter->first,iter->second, this->m_param)
	    );
	this->m_set.insert(make_pair(id1,id2));
	count_l++;
      }
    }
  }
  m_update_timestamp = this->m_ppa->getTimeStamp();

  console.Debug() << "added " << count_l << " pairs to ParallelInteractionStorage_ED\n";
  //  std::cout << "end ParallelInteractionStorage_ED::update at node " << m_comm.rank() << std::endl << std::flush;

	return res;
}


template<typename T,typename InteractionType>
void ParallelInteractionStorage_ED<T,InteractionType>::calcHeatFrict()
{
  console.Debug()
    << "calculating "
    << this->m_interactions.size()
    << " frictional heatings\n" ;

  for(
    typename list<InteractionType>::iterator it = this->m_interactions.begin();
    it != this->m_interactions.end();
    it++
  ){
    it->calcHeatFrict();
  }
}

template<typename P,typename InteractionType>
void ParallelInteractionStorage_ED<P,InteractionType>::calcHeatTrans()
{
  console.Debug()
    << "calculating "
    << this->m_interactions.size()
    << " heat transfers\n" ;

  for(
    typename list<InteractionType>::iterator it = this->m_interactions.begin();
    it != this->m_interactions.end();
    it++
  ){
    it->calcHeatTrans();
  }
}

/*! 
   save checkpoint (i.e. restart) data
*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_ED<P,InteractionType>::saveCheckPointData(std::ostream &oStream)
{
  const std::string delim = "\n";

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
void ParallelInteractionStorage_ED<P,InteractionType>::loadCheckPointData(std::istream &iStream)
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
  
    // -- read bonds
    for(int i=0;i<nconn;i++){
      InteractionType new_bond;
      // read a bond
      new_bond.loadRestartData(iStream);
      // insert it into interaction storage
      this->tryInsert(new_bond);
    }
  }
}
