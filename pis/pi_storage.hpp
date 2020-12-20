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

#include "Foundation/vec3.h"
#include "Fields/ScalarInteractionFieldSlave.h"
#include "Fields/CheckedScalarInteractionFieldSlave.h"
#include "Fields/ScalarInteractionFieldSlaveTagged.h"
#include "Fields/CheckedScalarInteractionFieldSlaveTagged.h"
#include "Fields/VectorInteractionFieldSlave.h"

#include "ppa/src/pp_array.h"

template <typename I>
bool TParallelInteractionStorage<I>::InteractionIterator::isInner(const Iterator &it)
{
  return m_ppa->isInInner(it->getPosFirst());
}

template <typename I>
TParallelInteractionStorage<I>::InteractionIterator::InteractionIterator(
  Iterator begin,
  Iterator end,
  AParallelParticleArray* ppa
)
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

template <typename I>
bool TParallelInteractionStorage<I>::InteractionIterator::hasNext()
{
  return (m_numRemaining > 0);
}

template <typename I>
typename TParallelInteractionStorage<I>::InteractionIterator::Interaction &
TParallelInteractionStorage<I>::InteractionIterator::next()
{
  while (!isInner(m_it)) {
    m_it++;
  }
  Interaction &i = *m_it;
  m_it++;
  m_numRemaining--;
  return i;
}

template <typename I>
int TParallelInteractionStorage<I>::InteractionIterator::getNumRemaining()
{
  return m_numRemaining;
}

template <typename I>
typename TParallelInteractionStorage<I>::InteractionIterator
TParallelInteractionStorage<I>::getInnerInteractionIterator()
{
  return
    InteractionIterator(m_interactions.begin(), m_interactions.end(), m_ppa);
}


/*!
  For all interactions with the lower particle in the inner area of the ntable call a 
  function reading a value and return the results in a vector of <position,vaule> pairs

  \param rdf the function
*/
template <typename I> 
template <typename P> 
vector<pair<Vec3,P> > TParallelInteractionStorage<I>::forAllInnerInteractionsGetWithPos(P (I::*rdf)() const)
{
  vector<pair<Vec3,P> > res;

  for(typename list<I>::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++){
    Vec3 pos=iter->getPosFirst();
    if(m_ppa->isInInner(pos)) res.push_back(make_pair(iter->getPos(),((*iter).*rdf)())); 
  }

  return res;
}

/*!
  For all interactions with the lower particle in the inner area of the ntable call a 
  function reading a value and return the results in a vector of <<pos1,radius1,pos2,radius2,ipos>,value> groups

  \param rdf the function
*/
template <typename I> 
template <typename P> 
vector<pair<typename TParallelInteractionStorage<I>::Raw2Data,P> >
TParallelInteractionStorage<I>::forAllInnerInteractionsGetRaw2(P (I::*rdf)() const)
{
  vector<pair<Raw2Data,P> > res;

  for(typename list<I>::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++){
    if(m_ppa->isInInner(iter->getPosFirst())) {      
      const Raw2Data data = iter->getRaw2Data();
      res.push_back(pair<Raw2Data,P>(data,((*iter).*rdf)())); 
    }
  }

  return res;
}

/*!
  For all interactions with the lower particle in the inner area of the ntable call a 
  function reading a value and return the results in a vector of  <<pid1,pid2,pos1,pos2,ipos>,value> groups

  \param rdf the function
*/
template <typename I> 
template <typename P> 
vector<pair<typename TParallelInteractionStorage<I>::DataWithPosID,P> > 
TParallelInteractionStorage<I>::forAllInnerInteractionsGetDataWithPosID(P (I::*rdf)() const)
{
  vector<pair<DataWithPosID,P> > res;

  for(typename list<I>::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++){
    if(m_ppa->isInInner(iter->getPosFirst())) {      
      vector<int> ids=iter->getAllID();
      int id1=ids[0];
      int id2;
      if(ids.size()>=2) {
	id2=ids[1];
      } else {
	id2=-1;
      }
      const Raw2Data data = iter->getRaw2Data();
      Vec3 pos1=data.get<0>();
      Vec3 pos2=data.get<2>();
      Vec3 ipos=data.get<4>();
      res.push_back(pair<DataWithPosID,P>(DataWithPosID(id1,id2,pos1,pos2,ipos),((*iter).*rdf)())); 
    }
  }
  
  return res;
}

/*!
  For all interactions with the lower particle in the inner area of the ntable call a 
  function reading a value and return the results in a vector of <<ipos,pid1,pid2>,value> groups

  \param rdf the function
*/
template <typename I> 
template <typename P> 
vector<pair<typename TParallelInteractionStorage<I>::DataWithID,P> > 
TParallelInteractionStorage<I>::forAllInnerInteractionsGetDataWithID(P (I::*rdf)() const)
{
  vector<pair<DataWithID,P> > res;

  for(typename list<I>::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++){
    if(m_ppa->isInInner(iter->getPosFirst())) {      
      vector<int> ids=iter->getAllID();
      int id1=ids[0];
      int id2;
      if(ids.size()>=2) {
	id2=ids[1];
      } else {
	id2=-1;
      }
      Vec3 pos=iter->getPos();
      res.push_back(pair<DataWithID,P>(DataWithID(id1,id2,pos),((*iter).*rdf)())); 
    }
  }
  
  return res;
}

/*!
  For all interactions with the lower particle in the inner area of the ntable call a 
  function reading a value and return the results in a container
  particle ids

  \cont the container
  \param rdf the function
*/
template <typename I> 
template <typename P> 
void TParallelInteractionStorage<I>::forAllInnerInteractionsGet(P& cont,typename P::value_type (I::*rdf)()const)
{
  for(typename list<I>::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++){
    Vec3 pos=iter->getPosFirst();
    if(m_ppa->isInInner(pos)) cont.push_back(((*iter).*rdf)());  
  }
}


/*!
  For all interactions with the lower particle in the inner area of the ntable and one of the
  particles having the specified tag call a  function reading a value and return the results 
  in a vector of <position,value> pairs

  \param rdf the function
  \param tag the tag
  \param mask the mask used in tag comparison
*/
template <typename I> 
template <typename P> 
vector<pair<Vec3,P> > TParallelInteractionStorage<I>::forAllTaggedInnerInteractionsGetWithPos(P (I::*rdf)() const,int tag,int mask)
{
  vector<pair<Vec3,P> > res;

  for(typename list<I>::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++){
    Vec3 pos=iter->getPosFirst();
    if(iter->hasTag(tag,mask)){
      if(m_ppa->isInInner(pos)) res.push_back(make_pair(iter->getPos(),((*iter).*rdf)())); 
    }
  }

  return res;
}

/*!
  For all interactions with the lower particle in the inner area of the ntable and one of the
  particles having the specified tag call a function reading a value and return the results in a container

  \cont the container
  \param rdf the function
  \param tag the tag
  \param mask the mask used in tag comparison
*/
template <typename I> 
template <typename P> 
void TParallelInteractionStorage<I>::forAllTaggedInnerInteractionsGet(P& cont,typename P::value_type (I::*rdf)()const,int tag,int mask)
{
  for(typename list<I>::iterator iter=m_interactions.begin();
      iter!=m_interactions.end();
      iter++){
    Vec3 pos=iter->getPosFirst();
    if(iter->hasTag(tag,mask)){
      if(m_ppa->isInInner(pos)) cont.push_back(((*iter).*rdf)());  
    }
  }
}

/*!
  For all interactions a member functions is called which sets a parameter of the interacion.
    
  \param pname the name of the parameter to be set 
  \param val the value
*/
template <typename I>
void TParallelInteractionStorage<I>::forAllInteractionsSet(const string& pname,double val)
{
	//void (I::*stf)() = 
	typename I::ScalarSetFunction stf=I::getScalarSetFunction(pname);
	// check if name was valid, i.e. got a setter function
	if(stf==NULL){
		console.Critical() << "NULL stf !!!\n";
	} else {
		console.XDebug() << "got stf\n";
		for(typename list<I>::iterator iter=m_interactions.begin();
		  iter!=m_interactions.end();
		  iter++){
			((*iter).*stf)(val);  
		}
	}
}

/*!
  generate new scalar field saver from the PIS

  \param comm 
  \param fieldname
  \param is_checked 
  \param is_tagged
  \param tag
  \param mask
*/
template <typename I> 
AFieldSlave* TParallelInteractionStorage<I>::generateNewScalarFieldSlave(TML_Comm* comm,const string& fieldname,int is_checked,int is_tagged,int tag,int mask)
{
  InteractionFieldSlave<I>* new_ifs;
 

  if(is_checked==0){
    typename I::ScalarFieldFunction rdf=I::getScalarFieldFunction(fieldname);
    if(is_tagged==0){
      new_ifs=new ScalarInteractionFieldSlave<I>(comm,this,rdf);
    } else {
      new_ifs=new ScalarInteractionFieldSlaveTagged<I>(comm,this,rdf,tag,mask);
    }
  } else {
    typename I::CheckedScalarFieldFunction rdf=I::getCheckedScalarFieldFunction(fieldname);
    if(is_tagged==0){
      new_ifs=new CheckedScalarInteractionFieldSlave<I>(comm,this,rdf);
    } else {
      new_ifs=new CheckedScalarInteractionFieldSlaveTagged<I>(comm,this,rdf,tag,mask);
    }
  }

  return new_ifs;
}

/*!
  generate new vector field saver from the PIS

  \param comm 
  \param fieldname
  \param is_checked 
  \param is_tagged
  \param tag
  \param mask
*/
template <typename I> 
AFieldSlave* TParallelInteractionStorage<I>::generateNewVectorFieldSlave(TML_Comm* comm,const string& fieldname,int is_checked,int is_tagged,int tag,int mask)
{
  InteractionFieldSlave<I>* new_ifs = NULL;
 

  if(is_checked==0){
    typename I::VectorFieldFunction rdf=I::getVectorFieldFunction(fieldname);
    if(is_tagged==0){
      new_ifs=new VectorInteractionFieldSlave<I>(comm,this,rdf);
    } else {
      // new_ifs=new VectorInteractionFieldSlaveTagged<I>(comm,this,rdf,tag,mask);
    }
  } else {
    // typename I::CheckedVectorFieldFunction rdf=I::getCheckedVectorFieldFunction(fieldname);
    if(is_tagged==0){
      // new_ifs=new CheckedVectorInteractionFieldSlave<I>(comm,this,rdf);
    } else {
      // new_ifs=new CheckedVectorInteractionFieldSlaveTagged<I>(comm,this,rdf,tag,mask);
    }
  }

  return new_ifs;
}

