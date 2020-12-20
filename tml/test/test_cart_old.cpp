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
#include "cart_comm.h"

//--- System includes ---
#include <iostream>

//--- STL ---
#include <vector>

int main(int argc, char** argv)
{
  int data0,data1=0;
  std::vector<int> dims;
  std::vector<int> vdata0,vdata1;
  std::vector<bool> period;
  
  MPI_Init(&argc,&argv);
  
  TML_CommWorld worldcomm;
  
  dims.push_back(worldcomm.size());
  period.push_back(true);

  TML_CartComm cartcomm(&worldcomm,1,dims,period);

  data0=cartcomm.rank()+42;
  vdata0.push_back(cartcomm.rank()+69);
  vdata0.push_back(cartcomm.rank()+70);
  cartcomm.shift_cont(vdata0,vdata1,0,1,0);
  std::cout << "rank " << worldcomm.rank() << "  recveived cont.: \n";
  for(std::vector<int>::iterator iter=vdata1.begin();
      iter!=vdata1.end();
      iter++){
    std::cout << *iter << " ";
  }
  std::cout << std::endl;

  //cartcomm.shift(data0,data1,0,1,0);
  //std::cout << "rank " << worldcomm.rank() << " recveived: " << data1 << std::endl;
  
  MPI_Finalize();

  return 0;
}
