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

#include "test_pack.h"

//--- I/O includes ---
#include <iostream>
using std::cout;
using std::endl;
using std::flush;

//--- STL ---
#include <vector>
using std::vector;

//--- project includes ---
#include "packed_message.h"
#include "packed_multi_message.h"

// test 1 :
// pack and unpack some buildin data types (no communication)
bool test_pack(TML_Comm *comm, int rank)
{
  bool res=true;

  TML_Packed_Message msg(comm->comm(),4); // force grow()
  int j,i=42;
  double e,d=2.718;

  if(rank==0){
    // pack
    msg.pack(i);
    msg.pack(d);
    cout << "packed size: " << msg.size() << endl;
    // unpack
    msg.begin_unpack();
    msg.unpack(j);
    msg.unpack(e);
    cout << "unpacked : " << j << " , " << e << endl;
    // ckeck
    res=(j==42)&&(e=2.718);
  }
  return res;
}

// test 2
// pack/unpack some buildin data types (no communication), checked

// pack/unpack buildin data types to/from multi-message
bool test_pack_multi(TML_Comm *comm, int rank)
{
  bool res=true;

  TML_PackedMultiMessage msg(comm->comm(),4); // force grow()
  int j,i=42;
  double e,d=2.718;

  if(rank==0){
    // pack
    msg[1].pack(i);
    msg[1].pack(d);
    //cout << "packed size: " << msg.size() << endl;
    // unpack
    msg.begin_unpack(1);
    msg[1].unpack(j);
    msg[1].unpack(e);
    cout << "multi unpacked : " << j << " , " << e << endl;
    // ckeck
    res=(j==42)&&(e=2.718);
  }
  return res;
}

// test 3 :
// packed send and receive a container of builtin values
bool test_container_packed(TML_Comm *comm, int rank)
{  
  vector<int> vdata0;
  vector<int> vdata1;
  bool res=true;

  if(rank==0){
    vdata0.push_back(69);
    vdata0.push_back(70);
    comm->send_cont_packed(vdata0,1,false,3);
    
  }else if (rank==1){
    comm->receive_cont_packed(vdata1,0,false,3);
    std::cout << "rank " << rank << "  received cont.: \n";
    for(std::vector<int>::iterator iter=vdata1.begin();
	iter!=vdata1.end();
	iter++){
      std::cout << *iter << " ";
    }
    std::cout << std::endl;
    // check
    if(vdata1.size()==2){
      res=(vdata1[0]==69)&&(vdata1[1]==70);
    } else {
      res=false;
    }
  }

  return res;
}

// test 4:
// packed sendrecv container of builtin values
bool test_container_sendrecv_packed(TML_Comm *comm, int rank)
{ 
  vector<int> vdata0;
  vector<int> vdata1;
  bool res=true;

  if(rank==0){
    //vdata0.push_back(69);
    //vdata0.push_back(70);
    comm->sendrecv_cont_packed(vdata0,vdata1,MPI_PROC_NULL,1,0);
    cout << "rank " << rank << "  received packed cont.: \n";
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
    comm->sendrecv_cont_packed(vdata0,vdata1,0,MPI_PROC_NULL,0);
    cout << "rank " << rank << "  received packed cont.: \n";
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

// test 5:
// packed replace sendrecv container of builtin values
bool test_container_sendrecv_packed_replace(TML_Comm *comm, int rank)
{ 
  vector<int> vdata;
  bool res=true;

  if(rank==0){
    vdata.push_back(69);
    vdata.push_back(70);
    comm->sendrecv_cont_packed_replace(vdata,1,1,0);
    cout << "rank " << rank << "  received packed cont.: \n";
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
    comm->sendrecv_cont_packed_replace(vdata,0,0,0);
    cout << "rank " << rank << "  received packed cont.: \n";
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

// packed send/receive of user defined type 
// new pack/unpack functions
template<>
void TML_PackedMessageInterface::pack<pair<int,double> >(const pair<int,double>& p)
{
  append(p.first);
  append(p.second);
}

template<>
void TML_PackedMessageInterface::unpack<pair<int,double> >(pair<int,double>& p)
{
  p.first=pop_int();
  p.second=pop_double();
}

bool test_container_packed_user(TML_Comm *comm, int rank)
{  
  vector<pair<int,double> > vdata0;
  vector<pair<int,double> > vdata1;
  bool res=true;

  if(rank==0){
    vdata0.push_back(pair<int,double>(69,2.718));
    vdata0.push_back(pair<int,double>(42,3.14159));
    comm->send_cont_packed(vdata0,1,false,3);
    
  }else if (rank==1){
    comm->receive_cont_packed(vdata1,0,false,3);
    std::cout << "rank " << rank << "  received cont.: \n";
    for(std::vector<pair<int,double> >::iterator iter=vdata1.begin();
	iter!=vdata1.end();
	iter++){
      std::cout << "(" << iter->first << "," << iter->second << ") " ;
    }
    std::cout << std::endl;
    // check
    if(vdata1.size()==2){
      res=(vdata1[0].first==69)&&(vdata1[0].second==2.718)&&
	(vdata1[1].first==42)&&(vdata1[1].second==3.14159);
    } else {
      res=false;
    }
  }

  return res;
}

// do all tests in this group
bool test_group_pack(TML_Comm *comm, int rank)
{
  bool res=true;
  cout << "begin test group packed" << endl;
  // test 1
  if(test_pack(comm,rank)){
    cout << "pack/unpack test sucessfull" << endl;
  }else{
    res=false;
    cout << "pack/unpack test failed" << endl;
  }
  // multi message pack/unpack
  if(test_pack_multi(comm,rank)){
    cout << "multi message pack/unpack test sucessfull" << endl;
  }else{
    res=false;
    cout << "multi message pack/unpack test failed" << endl;
  }
  // test 3
  if(test_container_packed(comm,rank)){
    cout << "test_container_packed sucessfull" << endl;
  }else{
    res=false;
    cout << "test_container_packed  failed" << endl;
  }
  // test 3 user def
  if(test_container_packed_user(comm,rank)){
    cout << "test_container_packed_user sucessfull" << endl << flush;
  }else{
    res=false;
    cout << "test_container_packed_user failed" << endl<< flush;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  cout << "past barrier" << endl << flush;
  // test 4
  if(test_container_sendrecv_packed(comm,rank)){
    cout << "test_container_sendrecv_packed sucessfull" << endl << flush;
  }else{
    res=false;
    cout << "test_container_sendrecv_packed failed" << endl << flush;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // test 5
  if(test_container_sendrecv_packed_replace(comm,rank)){
    cout << "test_container_sendrecv_packed_replace sucessfull" << endl << flush;
  }else{
    res=false;
    cout << "test_container_sendrecv_packed_replace failed" << endl << flush;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  cout << "finished test group packed" << endl << flush;

  return res;
}
