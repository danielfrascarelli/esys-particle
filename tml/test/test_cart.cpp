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

#include "test_cart.h"

//--- TML ---
#include <cart_comm.h>

//--- I/O includes ---
#include <iostream>
using std::cout;
using std::endl;

//--- STL ---
#include <vector>
using std::vector;

// construct cartesian communicator
TML_CartComm cart_comm_build(TML_Comm *comm, bool circular)
{
  std::vector<int> dims;
  std::vector<bool> period;

  dims.push_back(comm->size());
  period.push_back(circular);

  return TML_CartComm(comm,1,dims,period);
}

// construct cartesian communicator without given dimensions
bool test_cart_comm_build(TML_Comm *comm)
{
  std::vector<bool> period;

  TML_CartComm *cc=new TML_CartComm(comm,3,vector<int>(3,0),vector<bool>(3,false));
  
  delete cc;

  return true;
}



// get coords
bool test_get_coords(TML_CartComm *comm,int rank)
{
  bool res=true;

  vector<int> mycoords=comm->get_coords(); // my own coords
  vector<int> othercoords=comm->get_coords(1); // coords of rank 1

  // check own coords
  if(mycoords.size()==1){
    cout << "own coords of rank " << rank << " : " << mycoords[0] << endl;
    if(mycoords[0]!=rank){
      res=false;
    }
  } else {
    res=false;
  }

  //check other coords
 if(othercoords.size()==1){
    cout << "coords of rank 1 from rank " << rank << " : " << othercoords[0] << endl;
    if(othercoords[0]!=1){
      res=false;
    }
  } else {
    res=false;
  }

  return res;
}

// test shift functions
// because they are based on the sendrecv* stuff, only one (shift_cont) is tested
bool test_shift(TML_CartComm *comm,int rank,bool circ)
{
  bool res=true;
  std::vector<int> vdata0,vdata1;

  // setup data
  vdata0.push_back(rank+69);
  vdata0.push_back(rank+70);

  // shift by 1 in direction 0
  comm->shift_cont(vdata0,vdata1,0,1,42);
  
  //check result
  if(circ){ // circular boudary cond.
    
  }else{ //open boundary cond.
    if(rank==0){ // shouldn't get any
      if(vdata1.size()!=0){
	res=false;
	cout << "rank 0 shouldn't have received anything" << endl;
      }  
    } else { // should get 69,70
      if(vdata1.size()==2){
	res=(vdata1[0]==69) && (vdata1[1]=70);
	cout << "rank 1 received : " << vdata1[0] << " " << vdata1[1] << endl;
      } else {
	cout << "rank 1 wrong recv. size " << vdata1.size() << endl;
	res=false;
      }
    }
  }
  
  return res;
}

// do all tests in this group
bool test_group_cart(TML_Comm *comm, int rank)
{
  bool res=true;
  
  // build 2 cartesian comms (circular/open)
  TML_CartComm circular_comm=cart_comm_build(comm,true);
  TML_CartComm open_comm=cart_comm_build(comm,false);
  
  // tests
  if(test_get_coords(&circular_comm,rank)){
    cout << "test_get_coords sucessfull" << endl;
  }else{
  res=false;
    cout << "test_get_coords failed" << endl;
  }

  if(test_cart_comm_build(comm)){
    cout << "test_cart_comm sucessfull" << endl;
  }else{
  res=false;
    cout << "test_cart_comm failed" << endl;
  }

  // shift with open bcond
  if(test_shift(&open_comm,rank,false)){
    cout << "test_shift with open boundarys sucessfull" << endl;
  }else{
  res=false;
    cout << "test_shift with open boundarys failed" << endl;
  }

  cout << "finished test group cartesian" << endl;

  return res;
  }
