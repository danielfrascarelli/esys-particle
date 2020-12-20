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

//--- MPI ---
#include <mpi.h>

//--- TML ---
#include "comm_world.h"

//--- Test groups ---
#include "test_comm.h"
#include "test_pack.h"
#include "test_cart.h"
#include "test_sc.h"

//--- System includes ---
#include <iostream>
using std::cout;
using std::endl;
using std::flush;

//--- STL ---
#include <vector>
using std::vector;



int main(int argc, char** argv)
{
  bool res=true;
  
  MPI_Init(&argc,&argv);
  
  TML_CommWorld worldcomm;

  int rank=worldcomm.rank();
  
  // normal communication
  if(test_group_comm(&worldcomm,rank)){
    cout << "TML_Comm tests sucessfull" << endl << flush;
  }else{
    res=false;
    cout << "TML_Comm tests failed" << endl << flush;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // packed communication
  if(worldcomm.size()==3){ // currently requires 3 nodes
    if(test_group_pack(&worldcomm,rank)){
      cout << "TML_Comm packed communication tests sucessfull" << endl << flush;
    }else{
      res=false;
      cout << "TML_Comm packed communication  tests failed" << endl << flush;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // cartesian communicator communication
  if(test_group_cart(&worldcomm,rank)){
    cout << "TML_CartComm communication tests sucessfull" << endl << flush;
  }else{
    res=false;
    cout << "TML_CartComm communication  tests failed" << endl << flush;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  // scatter/gather tests
  if(worldcomm.size()==3){ // currently requires 3 node
    if(test_group_sc(&worldcomm,rank)){
      cout << "scatter/gather tests sucessfull" << endl << flush;
    } else{
      res=false;
      cout << "scatter/gather tests failed" << endl << flush;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  

  MPI_Finalize();
  if(res){
    cout << rank << " finalized - all tests succeeded" << endl;
  } else {
    cout << rank << " finalized - some tests failed" << endl;
  }
  return 0;
}
