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
using std::make_pair;

// -- IO includes --
#include <iostream>

using std::cout;
using std::endl;
using std::cerr;

// --- project includes ---
#include "Foundation/triplet.h"
#include "Foundation/console.h"

using esys::lsm::triplet;

/*!
  constructor

  \param comm the TML communicator used for sending the data back to the master
  \param pis a pointer to the interaction storage
  \param rdf the particle member function to access the data
*/
template <typename T>
CheckedScalarInteractionFieldSlave<T>::CheckedScalarInteractionFieldSlave(TML_Comm* comm,TParallelInteractionStorage<T>* pis,typename T::CheckedScalarFieldFunction rdf):InteractionFieldSlave<T>(comm,pis)
{
  this->m_rdf=rdf;
}

/*!
  send full field data and position of the interaction

	\todo replace cerr by exception
*/
template <typename T>
void CheckedScalarInteractionFieldSlave<T>::SendDataFull()
{
  vector<pair<Vec3,pair<bool,double> > > raw_data; // position, [valid ?, data]
  vector<pair<Vec3,double> > data; // position, data

  // get raw data
  raw_data=this->m_pis->forAllInnerInteractionsGetWithPos(this->m_rdf);
  
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
  send full field data and position of the interaction

	\todo replace cerr by exception
*/
template <typename T>
void CheckedScalarInteractionFieldSlave<T>::SendDataFull2()
{
  vector<pair<esys::lsm::quintuple<Vec3,double,Vec3,double,Vec3>,pair<bool,double> > > raw_data;
  vector<pair<esys::lsm::quintuple<Vec3,double,Vec3,double,Vec3>,double> > data; 

  // get raw data
  raw_data=this->m_pis->forAllInnerInteractionsGetRaw2(this->m_rdf);
  
  // filter data
  for(vector<pair<esys::lsm::quintuple<Vec3,double,Vec3,double,Vec3>,pair<bool,double> > >::iterator iter=raw_data.begin();
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
void CheckedScalarInteractionFieldSlave<T>::SendDataSum()
{
  vector<pair<bool,double> > data_vec;

  // get data from interactions
  this->m_pis->forAllInnerInteractionsGet(data_vec,this->m_rdf);

  // sum data
  double sum=0.0;
  for(vector<pair<bool,double> >::iterator iter=data_vec.begin();
      iter!=data_vec.end();
      iter++){
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
void CheckedScalarInteractionFieldSlave<T>::SendDataMax()
{
  vector<pair<bool,double> > data_vec;

  // get data from interactions
  this->m_pis->forAllInnerInteractionsGet(data_vec,this->m_rdf);

  // sum data
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

/*!
  send data together with <id1,id2,pos> info
*/ 
template <typename T>
void CheckedScalarInteractionFieldSlave<T>::SendDataWithID()
{
  vector<pair<triplet<int,int,Vec3>, pair<bool,double> > > raw_data;
  vector<pair<triplet<int,int,Vec3>, double> > data;

  // get raw data
  raw_data=this->m_pis->forAllInnerInteractionsGetDataWithID(this->m_rdf);

  console.XDebug() << "got " << raw_data.size() << " raw data\n"; 

  // filter data
  for(vector<pair<triplet<int,int,Vec3>, pair<bool,double> > >::iterator iter=raw_data.begin();
      iter!=raw_data.end();
      iter++){
    if(iter->second.first){
      data.push_back(make_pair(iter->first,iter->second.second));
    }
  }

  console.XDebug() << "got " << data.size() << " filtered data\n";

  // send data to master
  this->m_comm->send_gather(data,0);
}
