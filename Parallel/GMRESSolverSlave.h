/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2014 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0          //
//                                                         //
/////////////////////////////////////////////////////////////


#ifndef __GMRESSOLVERSLAVE_H
#define __GMRESSOLVERSLAVE_H

//--- STL includes ---
#include <string>
using std::string;

// --- project includes ---
#include "Foundation/console.h"
#include "Foundation/vec3.h"

// --- TML includes ---
#include "tml/comm/comm.h"

class TML_Comm;

/*!
  \class GMRESSolverSlave
  \brief Slave class of GMRES Solver for solving sparse matrix of linear equations

  \author Qi Shao
  $Revision$
  $Date$
*/
class GMRESSolverSlave
{
 private:
  int nrowsI,nproc,id,nnz;
  int Nx,Ny,Nz,Number,my_nnz;
  vector<int> irow,icol,ia,ja;
  vector<double> var, A, B, b;

 protected:
  MPI_Comm *m_worker_comm;
  TML_Comm *m_tml_global_comm;

  void LocalMatrixSetting();
  int max_int_array(vector<int>, int);
  void localize(int ***, int **, int ***, int **);
  void gather_sync(int **, int *, int **, int *, double *);
  void sparse_matmult(double *, double *);
  void dotproduct(double *, double *, double *);

 public:
  void LocalMatrixSolving();

  GMRESSolverSlave(
    MPI_Comm*,
    TML_Comm*,
    vector<double>,
    vector<int>,
    vector<int>,
    vector<double>,
    int,
    int
  );
   ~GMRESSolverSlave();
};

#endif //__GMRESSOLVERSLAVE_H



