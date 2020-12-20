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

#include "test_sc.h"

//--- I/O includes ---
#include <iostream>
using std::cout;
using std::endl;
using std::flush;

//--- STL ---
#include <vector>
#include <set>
#include <map>

using std::vector;
using std::set;
using std::multimap;

// test scatter
bool test_scatter(TML_Comm *comm, int rank)
{ 
  bool res=true;
  if(rank==0){ // rank 0 is sender
    // prepare send data
    multimap<int,double> sdata;
    sdata.insert(pair<int,double>(1,3.5));
    sdata.insert(pair<int,double>(1,4.5));
    sdata.insert(pair<int,double>(2,5.5));
    sdata.insert(pair<int,double>(2,6.5));
    sdata.insert(pair<int,double>(2,7.5));
    comm->scatter(sdata);
  } else { // other ranks receive
    vector<double> rdata;
    comm->recv_scatter(rdata,0);
    // test
    if(rank==1){
      res=(rdata.size()==2)&&(rdata[0]==3.5)&&(rdata[1]==4.5);
    } else if(rank==2){
      res=(rdata.size()==3)&&(rdata[0]==5.5)&&(rdata[1]==6.5)&&(rdata[2]==7.5);
    }
  }
  return res;
}

// test gather
bool test_gather(TML_Comm *comm, int rank)
{ 
  vector<double> sdata;
  multimap<int,double> rdata;
  bool res=true;
  // setup data
  if(rank==0){
    sdata.push_back(2.2);
    sdata.push_back(3.3);
    comm->send_gather(sdata,2);
  } else if(rank==1){
    sdata.push_back(4.4);
    sdata.push_back(5.5);
    sdata.push_back(6.6);
    comm->send_gather(sdata,2);
  } else { // rank 2 == root -> receive
    comm->gather(rdata);
    // test
    res=(rdata.count(0)==2) && (rdata.count(1)==3) && (rdata.count(2)==0);
    if(res){
      multimap<int,double>::const_iterator iter;
      iter=rdata.find(0);
      res=res&&(iter->second==2.2);
      iter++;
      res=res&&(iter->second==3.3);
      iter=rdata.find(1);
      res=res&&(iter->second==4.4);
      iter++;
      res=res&&(iter->second==5.5);
      iter++;
      res=res&&(iter->second==6.6);
    } 
  }
  return res;
}

// test scatter packed
bool test_scatter_packed(TML_Comm *comm, int rank)
{ 
  bool res=true;
  if(rank==0){ // rank 0 is sender
    // prepare send data
    multimap<int,double> sdata;
    sdata.insert(pair<int,double>(1,3.5));
    sdata.insert(pair<int,double>(1,4.5));
    sdata.insert(pair<int,double>(2,5.5));
    sdata.insert(pair<int,double>(2,6.5));
    sdata.insert(pair<int,double>(2,7.5));
    comm->scatter_packed(sdata);
  } else { // other ranks receive
    vector<double> rdata;
    comm->recv_scatter_packed(rdata,0);
    // test
    if(rank==1){
      res=(rdata.size()==2)&&(rdata[0]==3.5)&&(rdata[1]==4.5);
    } else if(rank==2){
      res=(rdata.size()==3)&&(rdata[0]==5.5)&&(rdata[1]==6.5)&&(rdata[2]==7.5);
    }
  }
  return res;
}

// test broadcast
bool test_broadcast(TML_Comm *comm, int rank)
{
  bool res=true;

  if(rank==0){ // rank 0 is sender
    comm->broadcast(42); 
  } else {
    int data;
    comm->recv_broadcast(data,0); 
    // test
    res=(data==42);
  }
  return res;
}

// test broadcast_cont
bool test_broadcast_cont(TML_Comm *comm, int rank)
{
  bool res=true;
  vector<double> sdata;
  vector<double> rdata;

  if(rank==0){ // rank 0 is sender
    sdata.push_back(3.14159);
    sdata.push_back(1.41421);
    comm->broadcast_cont(sdata);
  } else {
    comm->recv_broadcast_cont(rdata,0);
    //test
    if(rdata.size()==2){
      res=(rdata[0]==3.14159) && (rdata[1]==1.41421);
    }else{
      res=false;
      cout << "process " << rank << "received wrong size " << rdata.size() << endl;
    }
  }
  return res;
}

// test broadcast_cont_packed
bool test_broadcast_cont_packed(TML_Comm *comm, int rank)
{
  bool res=true;
  vector<double> sdata;
  vector<double> rdata;

  if(rank==0){ // rank 0 is sender
    sdata.push_back(3.14159);
    sdata.push_back(1.41421);
    comm->broadcast_cont_packed(sdata);
  } else {
    comm->recv_broadcast_cont_packed(rdata,0);
    //test
    if(rdata.size()==2){
      res=(rdata[0]==3.14159) && (rdata[1]==1.41421);
    }else{
      res=false;
      cout << "process " << rank << "received wrong size " << rdata.size() << endl;
    }
  }
  return res;
}


bool test_group_sc(TML_Comm *comm, int rank)
{
  bool res=true;
  
  if(test_scatter(comm,rank)){
    cout << "test_scatter sucessfull" << endl;
  }else{
    res=false;
    cout << "test_scatter failed" << endl;
  }

  if(test_gather(comm,rank)){
    cout << "test_gather sucessfull" << endl;
  }else{
    res=false;
    cout << "test_gather failed" << endl;
  }

  if(test_scatter_packed(comm,rank)){
    cout << "test_scatter_packed sucessfull" << endl;
  }else{
    res=false;
    cout << "test_scatter_packed failed" << endl;
  }

  if(test_broadcast(comm,rank)){
    cout << "test_broadcast sucessfull" << endl;
  }else{
    res=false;
    cout << "test_broadcast failed" << endl;
  }
  
  if(test_broadcast_cont(comm,rank)){
    cout << "test_broadcast_cont sucessfull" << endl;
  } else{
    res=false;
    cout << "test_broadcast_cont failed" << endl;
  }
  
  if(test_broadcast_cont_packed(comm,rank)){
    cout << "test_broadcast_cont_packed sucessfull" << endl;
  } else{
    res=false;
    cout << "test_broadcast_cont_packed failed" << endl;
  }

  cout << "finished test_group_sc"<< endl;
  return res;
}
