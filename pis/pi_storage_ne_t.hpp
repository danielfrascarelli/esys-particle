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
ParallelInteractionStorage_NE_T<P,I>::ParallelInteractionStorage_NE_T(AParallelParticleArray* ppa,const typename I::ParameterType& param,int tag1, int mask1, int tag2, int mask2):ParallelInteractionStorage_NE<P,I>(ppa,param)
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
bool ParallelInteractionStorage_NE_T<T,InteractionType>::update()
{
  console.XDebug() << "ParallelInteractionStorage_NE_T::Update\n";
  int count_l=0;
  bool res=true;

  if(this->m_update_timestamp != this->m_ppa->getTimeStamp()){// m_ppa rebuild since last update 
    // clean out old interactions
    this->m_interactions.clear();
    this->m_set.erase(this->m_set.begin(),this->m_set.end());
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
	tv.push_back(iter->first->getID());
	tv.push_back(iter->second->getID());
	if(!this->isExcluded(tv)){
	  this->m_interactions.push_back(InteractionType(iter->first,iter->second,this->m_param));
	  this->m_set.insert(pair<int,int>(iter->first->getID(),iter->second->getID()));
	  count_l++; 
	}
      }
    }
  } else { // m_ppa not rebuild since last update -> just get additional interactions
    // get list  of pairs from m_ppa
    typename ParallelParticleArray<T>::PairListHandle plh =
      ((ParallelParticleArray<T>*)this->m_ppa)->getNewPairList();
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
	tv.push_back(iter->first->getID());
	tv.push_back(iter->second->getID());
	if(!this->isExcluded(tv)){
	  this->m_interactions.push_back(InteractionType(iter->first,iter->second,this->m_param));
	  this->m_set.insert(pair<int,int>(iter->first->getID(),iter->second->getID()));
	  count_l++;
	}
      }
    }
  }
  this->m_update_timestamp = this->m_ppa->getTimeStamp();

  console.XDebug() << "added " << count_l << " pairs to EIG\n";
  console.XDebug() << "end ParallelInteractionStorage_NE_T::Update\n";

  return res;
}
