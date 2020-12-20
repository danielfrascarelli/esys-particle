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
ParallelInteractionStorage_NE<P,I>::ParallelInteractionStorage_NE(AParallelParticleArray* ppa,const typename I::ParameterType& param):TParallelInteractionStorage<I>(ppa)
{
  m_param=param;
  m_update_timestamp=0;
}

template<typename P,typename InteractionType>
void  ParallelInteractionStorage_NE<P,InteractionType>::addExIG(AParallelInteractionStorage* eg)
{
  m_exIG.push_back(eg); 
}

/*!
  check if a potential interaction is in one of the excluding IGs

  \param tv the potential interaction
*/
template<typename T,typename InteractionType>
bool ParallelInteractionStorage_NE<T,InteractionType>:: isExcluded(const vector<int> tv)
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

/*!
  Update interactions. Do full dynamic search.
*/
template<typename T,typename InteractionType>
bool ParallelInteractionStorage_NE<T,InteractionType>::update()
{
  console.XDebug() << "ParallelInteractionStorage_NE::Update\n";
  int count_l=0;
  bool res=true;

  if(m_update_timestamp != this->m_ppa->getTimeStamp()){// m_ppa rebuild since last update 
    // clean out old interactions
    this->m_interactions.clear();
    m_set.erase(m_set.begin(),m_set.end());
    // get list  of pairs from m_ppa
    typename ParallelParticleArray<T>::PairListHandle plh =
      ((ParallelParticleArray<T>*)this->m_ppa)->getFullPairList();
    // generate interactions from pairs
    for(typename ParallelParticleArray<T>::PairListIterator iter=plh->begin();
        iter!=plh->end();
        iter++){
      // check vs. ExIG
      vector<int> tv;
      tv.push_back(iter->first->getID());
      tv.push_back(iter->second->getID());
      if(!isExcluded(tv)){
	this->m_interactions.push_back(InteractionType(iter->first,iter->second,m_param));
	m_set.insert(pair<int,int>(iter->first->getID(),iter->second->getID()));
	count_l++;
      }
    }
  } else { // m_ppa not rebuild since last update -> just get additional interactions
    // get list  of pairs from m_ppa
    typename ParallelParticleArray<T>::PairListHandle plh =
      ((ParallelParticleArray<T>*)this->m_ppa)->getNewPairList();
    //cout << "got NewPairList: ";
    for(typename ParallelParticleArray<T>::PairListIterator iter=plh->begin();
        iter!=plh->end();
        iter++){
      // cout << iter->first->getID() << "-" << iter->second->getID() << endl;
      // check vs. ExIG
      vector<int> tv;
      tv.push_back(iter->first->getID());
      tv.push_back(iter->second->getID());
      if(!isExcluded(tv)){
	this->m_interactions.push_back(InteractionType(iter->first,iter->second,m_param));
	m_set.insert(pair<int,int>(iter->first->getID(),iter->second->getID()));
	count_l++;
      }
    }
  }
  m_update_timestamp = this->m_ppa->getTimeStamp();

  console.XDebug() << "added " << count_l << " pairs to EIG\n";
  console.XDebug() << "end ParallelInteractionStorage_NE::Update\n";

	return res;
}


/*!
  \warning evil hack, only checks 1st & 2nd id -> change from pair<int,int> to vector<int>
*/
template<typename P,typename InteractionType>
bool ParallelInteractionStorage_NE<P,InteractionType>::isIn(const vector<int>& pids)
{
  bool res;

  res=m_set.find(make_pair(pids[0],pids[1]))!=m_set.end();

  return res;
}

/*!

*/
template<typename P,typename InteractionType>
void ParallelInteractionStorage_NE<P,InteractionType>::calcForces()
{
  console.Debug()
    << "calculating "
    << this->m_interactions.size()
    << " interaction forces\n";

  for(
    typename list<InteractionType>::iterator it = this->m_interactions.begin();
    it != this->m_interactions.end();
    it++
  ){
    it->calcForces();
  }
}

template<typename P,typename InteractionType>
void ParallelInteractionStorage_NE<P,InteractionType>::calcHeatTrans()
{
  console.Debug()
    << "calculating "
    << this->m_interactions.size()
    << " interaction heat transfers\n" ;

  for(
    typename list<InteractionType>::iterator it = this->m_interactions.begin();
    it != this->m_interactions.end();
    it++
  ){
    it->calcHeatTrans();
  }
}
