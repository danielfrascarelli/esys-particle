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

#include "test_comm.h"

//--- I/O includes ---
#include <iostream>
using std::cout;
using std::endl;
using std::flush;

//--- STL ---
#include <vector>
#include <set>
#include <utility>

using std::vector;
using std::set;
using std::pair;

//--- other includes ---
#include "vec3.h"
#include "vec3_mpi.h"

// test communicator construction
bool test_const(TML_Comm *comm, int rank)
{
  bool res=true;
  vector<int> ids;
  ids.push_back(0);

  TML_Comm newcomm=comm->include(ids);
  cout << "rank " << rank << "past 1" << endl << flush;
  int nrank=newcomm.rank();
  cout << "rank " << rank << "past 2" << endl << flush;

  if(rank==0){
    res=(nrank==0);
  } else {
    res=(nrank==MPI_UNDEFINED);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  return res;
}

// send and receive simple data types
bool test_simple(TML_Comm *comm, int rank)
{
  int data0,data1=0;
  bool res=true;

  if(rank==0){
    data0=69;
    comm->send(data0,1,0);
  } else if(rank==1){
    comm->receive(data1,0,0);
    cout << "rank " << rank << " received: " << data1 << endl;
    res=(data1==69);
  }

  return res;
}

// test 2 :
// sendrecv simple data
bool test_simple_sendrecv(TML_Comm *comm, int rank)
{
  int data0,data1=0;
  bool res=true;

  if(rank==0){
    data0=69;
    comm->sendrecv(data0,data1,1,1,2);
    cout << "rank " << rank << " received: " << data1 << endl;
    res=(data1==42);
    res=true;
  } else if(rank==1){
    data0=42;
    comm->sendrecv(data0,data1,0,0,2);
    cout << "rank " << rank << " received: " << data1 << endl;
    res=(data1==69);
  }

  return res;
}

// test 3 :
// send and receive array data
bool test_array(TML_Comm *comm, int rank)
{
  return true;
}

// test 4 :
// sendrecv array data
bool test_array_sendrecv(TML_Comm *comm, int rank)
{
  return true;
}



// test 5 :
// send and receive container
bool test_container(TML_Comm *comm, int rank)
{ 
  vector<int> vdata0;
  vector<int> vdata1;
  bool res=true;

  // test with sequence container (vector)
  if(rank==0){
    vdata0.push_back(69);
    vdata0.push_back(70);
    comm->send_cont(vdata0,1,0);    
  }else if (rank==1){
    vdata0.push_back(71);
    comm->receive_cont(vdata1,0,0);
    std::cout << "rank " << rank << "  received cont.: \n";
    for(std::vector<int>::iterator iter=vdata1.begin();
	iter!=vdata1.end();
	iter++){
      std::cout << *iter << " ";
    }
    std::cout << std::endl;
    // check
    if(vdata1.size()==3){
      res=(vdata1[0]==71)&&(vdata1[1]==69)&&(vdata1[2]==70);
    } else {
      res=false;
    }
  }

  // test with associative container (set)  
  set<double> mdata0;
  set<double> mdata1;
  if(rank==0){
    mdata0.insert(69.0);
    mdata0.insert(70.0);
    comm->send_cont(mdata0,1,0);    
  }else if (rank==1){
    comm->receive_cont(mdata1,0,0);
    std::cout << "rank " << rank << "  received cont.: \n";
    for(std::set<double>::iterator iter=mdata1.begin();
	iter!=mdata1.end();
	iter++){
      std::cout << *iter << " ";
    }
    std::cout << std::endl;
    // check
    if(mdata1.size()==2){
      res=(mdata1.find(69)!=mdata1.end())&&(mdata1.find(70)!=mdata1.end());
    } else {
      res=false;
    }
  }
  return res;
}

// test 6 :
// sendrecv container
bool test_container_sendrecv (TML_Comm *comm, int rank)
{ 
  vector<int> vdata0;
  vector<int> vdata1;
  bool res=true;

  if(rank==0){
    vdata0.push_back(69);
    vdata0.push_back(70);
    comm->sendrecv_cont(vdata0,vdata1,1,1,0);
    cout << "rank " << rank << "  received cont.: \n";
    for(std::vector<int>::iterator iter=vdata1.begin();
	iter!=vdata1.end();
	iter++){
      cout << *iter << " ";
    }
    // check 
    if(vdata1.size()==3){
      res=(vdata1[0]==42)&&(vdata1[1]==43)&&(vdata1[2]==44);
    } else {
      res=false;
    }
  }else if (rank==1){
    vdata0.push_back(42);
    vdata0.push_back(43);
    vdata0.push_back(44);
    comm->sendrecv_cont(vdata0,vdata1,0,0,0);
    cout << "rank " << rank << "  received cont.: \n";
    for(std::vector<int>::iterator iter=vdata1.begin();
	iter!=vdata1.end();
	iter++){
      cout << *iter << " ";
    }
    cout << endl;
    // check 
    if(vdata1.size()==2){
      res=(vdata1[0]==69)&&(vdata1[1]==70);
    } else {
      res=false;
    }
  }

  return res;
}

// test 7 :
// sendrecv_replace container
bool test_container_sendrecv_replace (TML_Comm *comm, int rank)
{ 
  vector<int> vdata;
  bool res=true;

  if(rank==0){
    vdata.push_back(69);
    vdata.push_back(70);
    comm->sendrecv_cont_replace(vdata,1,1,0);
    cout << "rank " << rank << "  received cont.: \n";
    for(std::vector<int>::iterator iter=vdata.begin();
	iter!=vdata.end();
	iter++){
      cout << *iter << " ";
    }
    // check 
    if(vdata.size()==3){
      res=(vdata[0]==42)&&(vdata[1]==43)&&(vdata[2]==44);
    } else {
      res=false;
    }
  }else if (rank==1){
    vdata.push_back(42);
    vdata.push_back(43);
    vdata.push_back(44);
    comm->sendrecv_cont_replace(vdata,0,0,0);
    cout << "rank " << rank << "  received cont.: \n";
    for(std::vector<int>::iterator iter=vdata.begin();
	iter!=vdata.end();
	iter++){
      cout << *iter << " ";
    }
    cout << endl;
    // check 
    if(vdata.size()==2){
      res=(vdata[0]==69)&&(vdata[1]==70);
    } else {
      res=false;
    }
  }

  return res;
}


// send and receive container of complex data (Vec3)
bool test_container_vec3(TML_Comm *comm, int rank)
{ 
  vector<Vec3> vdata0;
  vector<Vec3> vdata1;
  bool res=true;

  // test with sequence container (vector)
  if(rank==0){
    vdata0.push_back(Vec3(1.0,2.0,3.0));
    vdata0.push_back(Vec3(4.0,5.0,6.0));
    comm->send_cont(vdata0,1,0);    
  }else if (rank==1){
    comm->receive_cont(vdata1,0,0);
    std::cout << "rank " << rank << "  received cont.: \n";
    for(std::vector<Vec3>::iterator iter=vdata1.begin();
	iter!=vdata1.end();
	iter++){
      std::cout << *iter << " ";
    }
    std::cout << std::endl;
    // check
    if(vdata1.size()==2){
      res=(vdata1[0].X()=1.0)&&(vdata1[0].Y()==2.0)&&(vdata1[0].Z()==3.0)&&
	(vdata1[1].X()=4.0)&&(vdata1[1].Y()==5.0)&&(vdata1[1].Z()==6.0);
    } else {
      res=false;
    }
  }


  return res;
}

// send and receive pairs
bool test_pair(TML_Comm *comm, int rank)
{
  pair<int,double> data0,data1;
  bool res=true;

  if(rank==0){
    data0=make_pair(42,2.71828);
    comm->send(data0,1,0);
  } else if(rank==1){
    comm->receive(data1,0,0);
    cout << "rank " << rank << " received: [" << data1.first 
	 << "," << data1.second <<"]"<< endl;
    res=(data1.first==42)&&(data1.second==2.71828);
    cout << "diff: " << data1.second-2.71828 << endl;
  }

  return res;
}

// send and receive map
bool test_container_map(TML_Comm *comm, int rank)
{ 
  map<int,double> vdata0;
  map<int,double> vdata1;
  bool res=true;

  // test with sequence container (vector)
  if(rank==0){
    vdata0[2]=2.71828;
    vdata0[5]=3.14159;
    comm->send_cont(vdata0,1,0);    
  }else if (rank==1){
    comm->receive_cont(vdata1,0,0);
    cout << "rank " << rank << "  received map.: \n";
    for(map<int,double>::iterator iter=vdata1.begin();
	iter!=vdata1.end();
	iter++){
      cout << "[" << iter->first << ","  << iter->second << "]  ";
    }
    cout << endl;
    // check
    if(vdata1.size()==2){
      res=(vdata1[2]==2.71828)&&(vdata1[5]==3.14159);
    } else {
      res=false;
    }
  }
  return res;
}

// test summation and distribution to all
bool test_sum_all(TML_Comm *comm, int rank)
{ 
  bool res=true;
  double data=rank;
  double sum;
  
  sum=comm->sum_all(data);
  
  // check
  int n=comm->size();
  double expected_sum=double((n*(n-1))/2);
  cout << "comm size= " << n << " expected sum = " << expected_sum << " sum= " << sum << endl;
  res=(sum==expected_sum);

  return res;
}

bool test_group_comm(TML_Comm *comm, int rank)
{
  bool res=true;

  cout << "begin test group comm" << endl;

  if(test_const(comm,rank)){
    cout << "test_const sucessfull" << endl;
  }else{
    res=false;
    cout << "test_const failed" << endl;
  }

  // test 1
  if(test_simple(comm,rank)){
    cout << "test_simple sucessfull" << endl;
  }else{
    res=false;
    cout << "test_simple failed" << endl;
  }
 
 // test pair
  if(test_pair(comm,rank)){
    cout << "test_pair sucessfull" << endl;
  }else{
    res=false;
    cout << "test_pair failed" << endl;
  }

  // test 2
  if(test_simple_sendrecv(comm,rank)){
    cout << "test_simple_sendrecv sucessfull" << endl;
  }else{
    res=false;
    cout << "test_simple_sendrecv failed" << endl;
  }

  // test 5
  if(test_container(comm,rank)){
    cout << "test_container sucessfull" << endl;
  }else{
    res=false;
    cout << "test_container failed" << endl;
  }
  
  if(test_container_vec3(comm,rank)){
    cout << "test_container_vec3 sucessfull" << endl;
  }else{
    res=false;
    cout << "test_container_vec3 failed" << endl;
  }
  
  if(test_container_map(comm,rank)){
    cout << "test_container_map sucessfull" << endl;
  }else{
    res=false;
    cout << "test_container_map failed" << endl;
  }
  
  // test 6
  if(test_container_sendrecv(comm,rank)){
    cout << "test_container_sendrecv sucessfull" << endl;
  }else{
    res=false;
    cout << "test_container_sendrecv failed" << endl;
  }
  
  // test 7
  if(test_container_sendrecv_replace(comm,rank)){
    cout << "test_container_sendrecv_replace sucessfull" << endl;
  }else{
    res=false;
    cout << "test_container_sendrecv_replace failed" << endl;
  }

  // test sum_all
  if(test_sum_all(comm,rank)){
    cout << "test_sum_all sucessfull" << endl;
  }else{
    res=false;
    cout << "test_sum_all failed" << endl;
  }
 
  cout << "finished test group comm" << endl;

  return res;
}
