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

#include "ppa/src/pp_array.h"
//----------------------------------------------
//       CElasticInteractionGroup functions
//----------------------------------------------
template<class T>
CRotThermElasticInteractionGroup<T>::CRotThermElasticInteractionGroup()
{
  m_exIG=NULL;
  this->m_update_timestamp=0;
}

template<class T>
CRotThermElasticInteractionGroup<T>::CRotThermElasticInteractionGroup(const CRotThermElasticIGP* Param)
{
  m_exIG=NULL;
  this->m_kr=Param->getSpringConst(); 
  m_diffusivity = Param->getDiffusivity()  ;
  this->m_update_timestamp=0;
}

template<class T>
void CRotThermElasticInteractionGroup<T>::setParam(const CRotThermElasticIGP* Param)
{
  this->m_kr=Param->getSpringConst();
  m_diffusivity = Param->getDiffusivity()  ;
}

template<class T>
void CRotThermElasticInteractionGroup<T>::calcForces()
{

  console.Debug() << "calculating " << m_interactions.size() << " elastic forces\n" ;

  for(vector<CRotThermElasticInteraction>::iterator it=m_interactions.begin();it!=m_interactions.end();it++){
    it->calcForces();
  }
}

/*!
  Update the interactions from the neighbortable. If the ParallelParticleArray has been rebuild since the last update, the interactiongroup is cleared and rebuild from scratch. Otherwise only new interactions are added.
*/
template<class T>
void CRotThermElasticInteractionGroup<T>::Update(ParallelParticleArray<T>* PPA)
{
  console.XDebug() << "CElasticInteractionGroup::Update\n";
  int count_l=0;

  if(this->m_update_timestamp!=PPA->getTimeStamp()){// PPA rebuild since last update 
    // clean out old interactions
    m_interactions.erase(m_interactions.begin(),m_interactions.end());
    this->m_set.erase(this->m_set.begin(),this->m_set.end());
    // get list  of pairs from PPA
    typename ParallelParticleArray<T>::PairListHandle plh=PPA->getFullPairList();
    // generate interactions from pairs
    for(typename ParallelParticleArray<T>::PairListIterator iter=plh->begin();
	iter!=plh->end();
	iter++){
      // check vs. ExIG
      vector<int> tv;
      tv.push_back(iter->first->getID());
      tv.push_back(iter->second->getID());
      if(m_exIG!=NULL){
	if(!this->m_exIG->isIn(tv)){
	  m_interactions.push_back(CRotThermElasticInteraction(iter->first,iter->second,m_k));
	  this->m_set.insert(pair<int,int>(iter->first->getID(),iter->second->getID()));
	  console.XDebug()<<"adding pair: " << iter->first->getID() << " - " 
			  << iter->second->getID() << "\n";
	  count_l++;
	} else {
	  console.XDebug()<<"not adding pair: " << iter->first->getID() << " - " 
			  << iter->second->getID() << "\n";
	}
      } else {
	m_interactions.push_back(CRotThermElasticInteraction(iter->first,iter->second,m_k));
	this->m_set.insert(pair<int,int>(iter->first->getID(),iter->second->getID()));
	console.XDebug()<<"adding pair: " << iter->first->getID() << " - " 
			  << iter->second->getID() << "\n";
      }
    }
  } else { // PPA not rebuild since last update -> just get additional interactions
    // get list  of pairs from PPA
    typename ParallelParticleArray<T>::PairListHandle plh=PPA->getNewPairList();
    for(typename ParallelParticleArray<T>::PairListIterator iter=plh->begin();
	iter!=plh->end();
	iter++){
      // check vs. ExIG
      vector<int> tv;
      tv.push_back(iter->first->getID());
      tv.push_back(iter->second->getID());
      if(this->m_exIG!=NULL){
	if(!this->m_exIG->isIn(tv)){
	  m_interactions.push_back(CRotThermElasticInteraction(iter->first,iter->second,m_k));
	  this->m_set.insert(pair<int,int>(iter->first->getID(),iter->second->getID()));
	  console.XDebug()<<"adding pair: " << iter->first->getID() << " - " 
			  << iter->second->getID() << "\n";
	  count_l++;
	} else {
	  console.XDebug()<<"not adding pair: " << iter->first->getID() << " - " 
			  << iter->second->getID() << "\n";
	}
      } else {
	m_interactions.push_back(CRotThermElasticInteraction(iter->first,iter->second,m_k));
	this->m_set.insert(pair<int,int>(iter->first->getID(),iter->second->getID()));
	console.XDebug()<<"adding pair: " << iter->first->getID() << " - " 
			<< iter->second->getID() << "\n";
      }
    }
  }
  this->m_update_timestamp=PPA->getTimeStamp();
  console.XDebug() << "added " << count_l << " pairs to EIG\n";
 
  console.XDebug() << "end CElasticInteractionGroup::Update\n";
}

template<class T>
ostream& operator<<(ostream& ost,const CRotThermElasticInteractionGroup<T>& E)
{
  ost << "CElasticInteractionGroup : \n";
  for(vector<CRotThermElasticInteraction>::const_iterator it=E.m_interactions.begin();it!=E.m_interactions.end();it++){
    ost << *it << " , " ;
  }
  ost << "exchange list: \n";
  for(set<pair<int,int> >::const_iterator vit=E.m_exchg_list.begin();vit!=E.m_exchg_list.end();vit++){
    ost << "[ " << vit->first << " from " << vit->second << " ]  ,  ";
  }
  return ost;
}
