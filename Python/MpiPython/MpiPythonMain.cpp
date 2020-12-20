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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>
#include <Python.h>
#include <iostream>

#include "Foundation/PathUtil.h"
#include <stdexcept>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/path.hpp>
#include <patchlevel.h>

using namespace boost;


//--project includes--
#include "Foundation/console.h"
#include "Parallel/SubLatticeControler.h"


int main( int argc, char **argv ) 
{
  esys::lsm::setPathEnv(argc, argv);
  
    int status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) {
      std::cerr << argv[0] << ": MPI_Init failed, exiting." << std::endl;
      return status;
    } 
    // get rank
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    if(myrank==0){    // if rank==0 --> master
      // start python
      //     std::cout << "master start\n";
      Py_Initialize();
      
#if PY_VERSION_HEX >= 0x03000000
      wchar_t** wargv = new wchar_t*[argc+1];
      for (int i=0;i<argc;++i)
      {
        int len = strlen(argv[i]);
        wargv[i] = new wchar_t[len+1];
        for (int j=0; j<len;++j)
        {
          wargv[i][j] = wchar_t(argv[i][j]);
        }
        wargv[i][len] = 0;
      }
      wargv[argc] = 0;
      status = Py_Main(argc, wargv);
#else
      status = Py_Main(argc, argv);
#endif
      
      Py_Finalize();
    } else { // if rank!=0 --> slave
      // start worker
      //      std::cout << "slave start\n";

    CSubLatticeControler SLC;

    SLC.initMPI();
    SLC.run();
    }
    
    MPI_Finalize();
 
    return status;
}

