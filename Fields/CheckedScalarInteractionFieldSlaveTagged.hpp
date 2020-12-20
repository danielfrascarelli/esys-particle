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

//-- STL includes --
#include <vector>
#include <utility>

using std::vector;
using std::pair;

// -- IO includes --
#include <iostream>

using std::cout;
using std::endl;

/*!
  constructor
 
  \param comm the TML communicator used for sending the data back to the master
  \param pis a pointer to the interaction storage
  \param rdf the particle member function to access the data
*/
template <typename T>
CheckedScalarInteractionFieldSlaveTagged<T>::CheckedScalarInteractionFieldSlaveTagged(TML_Comm* comm,TParallelInteractionStorage<T>* pis,typename T::CheckedScalarFieldFunction rdf,int tag,int mask):CheckedScalarInteractionFieldSlave<T>(comm,pis,rdf)
{
	m_tag=tag;
	m_mask=mask;
}

/*!
  send full field data and position of the interaction
*/
template <typename T>
void CheckedScalarInteractionFieldSlaveTagged<T>::SendDataFull()
{
  vector<pair<Vec3,pair<bool,double> > > raw_data;
  vector<pair<Vec3,double> > data;

  // get raw field data from interaction storage
  raw_data=this->m_pis->forAllTaggedInnerInteractionsGetWithPos(this->m_rdf,m_tag,m_mask);

  // filter data
  for(vector<pair<Vec3,pair<bool,double> > >::iterator iter=raw_data.begin();
      iter!=raw_data.end();
      iter++){
    if(iter->second.first){
      data.push_back(make_pair(iter->first,iter->second.second));
    }
  }
  
  // send data to master
  this->m_comm->send_gather(data,0);
}

/*!
  send sum only
*/
template <typename T>
void CheckedScalarInteractionFieldSlaveTagged<T>::SendDataSum()
{
	vector<pair<bool,double> >data_vec;

	// get data from interactions
	this->m_pis->forAllTaggedInnerInteractionsGet(data_vec,this->m_rdf,m_tag,m_mask);

	// sum data
	double sum=0.0;
	for(vector<pair<bool,double> >::iterator iter=data_vec.begin();
	        iter!=data_vec.end();
	        iter++)
	{
		if(iter->first) sum+=iter->second;
	}

	vector<double> sum_vec;
	sum_vec.push_back(sum);
	this->m_comm->send_gather(sum_vec,0);
}

/*!
  send maximum only
*/
template <typename T>
void CheckedScalarInteractionFieldSlaveTagged<T>::SendDataMax()
{
	vector<pair<bool,double> >data_vec;

	// get data from interactions
	this->m_pis->forAllTaggedInnerInteractionsGet(data_vec,this->m_rdf,m_tag,m_mask);

	// get max 
	double max;
	bool is_set=false;
	for(vector<pair<bool,double> >::iterator iter=data_vec.begin();
	    iter!=data_vec.end();
	    iter++){
	 if(iter->first) {
	   if(is_set){
	   	 max=(iter->second > max) ? iter->second : max;
	  } else {
	     max=iter->second;
	     is_set=true;
	  	 }
		}
	}

	vector<double> max_vec;
	max_vec.push_back(max);
	this->m_comm->send_gather(max_vec,0);
}
