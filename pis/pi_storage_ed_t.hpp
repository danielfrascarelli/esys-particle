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
ParallelInteractionStorage_ED_T<P,I>::ParallelInteractionStorage_ED_T(AParallelParticleArray* ppa,const typename I::ParameterType& param,int tag1, int mask1, int tag2, int mask2):ParallelInteractionStorage_ED<P,I>(ppa,param)
{
  if(tag1<=tag2){ // sort tags so that m_tag1<=m_tag2
    m_tag1=tag1;
    m_mask1=mask1;
    m_tag2=tag2;
    m_mask2=mask2;
  } else {
    m_tag1=tag2;
    m_mask1=mask2;
    m_tag2=tag1;
    m_mask2=mask1;
  }
}


/*!
  Update interactions. Do full dynamic search.
*/
template<typename T,typename InteractionType>
bool ParallelInteractionStorage_ED_T<T,InteractionType>::update()
{
  console.XDebug() << "ParallelInteractionStorage_ED_T::update at node " << this->m_comm.rank() << "\n";

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
      //--- check particle tags ---
      // get tags
      int t1=iter->first->getTag();
      int t2=iter->second->getTag();
      // sort tags
      if(t1>t2){
	int th=t1;
	t1=t2;
	t2=th;
      }
      // tags fit -> go on
      if(((t1 & m_mask1)==(m_tag1 & m_mask1)) && ((t2 & m_mask2)==(m_tag2 & m_mask2))){
	// check vs. ExIG
	vector<int> tv;
	// ids in pair
	int id1=iter->first->getID();
	int id2=iter->second->getID();
	tv.push_back(id1);
	tv.push_back(id2);

	if((!(this->isExcluded(tv)))&&(!this->isIn(tv))){  // if not already in or in ExIG
	  this->m_interactions.push_back(InteractionType(iter->first,iter->second,this->m_param));
	  this->m_set.insert(make_pair(id1,id2));
	  count_l++;
	}
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
      //--- check particle tags ---
      // get tags
      int t1=iter->first->getTag();
      int t2=iter->second->getTag();
      // sort tags
      if(t1>t2){
	int th=t1;
	t1=t2;
	t2=th;
      }
      // tags fit -> go on
      if(((t1 & m_mask1)==(m_tag1 & m_mask1)) && ((t2 & m_mask2)==(m_tag2 & m_mask2))){      
	// check vs. ExIG
	vector<int> tv;
	// ids in pair
	int id1=iter->first->getID();
	int id2=iter->second->getID();
	tv.push_back(id1);
	tv.push_back(id2); 
	if((!(this->isExcluded(tv)))&&(!(this->isIn(tv)))) {
	  this->m_interactions.push_back(InteractionType(iter->first,iter->second, this->m_param));
	  this->m_set.insert(make_pair(id1,id2));
	  count_l++;
	}
      }
    }
  }
  this->m_update_timestamp = this->m_ppa->getTimeStamp();

  console.Debug() << "added " << count_l << " pairs to ParallelInteractionStorage_ED_T\n";

  return res;
}
