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


#ifndef __GMRESSOLVERMASTER_H
#define __GMRESSOLVERMASTER_H

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
  \class GMRESSolverMaster
  \brief Master class of GMRES Solver for solving sparse matrix of linear equations

  \author Qi Shao
  $Revision$
  $Date$
*/
class GMRESSolverMaster
{
 private:
  int Nx,Ny,Nz,Number,nnz;
  vector<int> irow,icol;
  vector<double> W, E, N, S, D, U, C, B;
  vector <double> var;
  vector<Vec3> index;
  int m_rank;

 protected:
  MPI_Comm *m_comm;
  TML_Comm *m_tml_comm;

  void GlobalMatrixSetting();

 public:
  vector<pair<Vec3,double> > MatrixSolving();

  GMRESSolverMaster(
    MPI_Comm*,
    TML_Comm*,
    vector<pair<Vec3,double> >,
    vector<pair<Vec3,double> >,
    vector<pair<Vec3,double> >,
    vector<pair<Vec3,double> >,
    vector<pair<Vec3,double> >,
    vector<pair<Vec3,double> >,
    vector<pair<Vec3,double> >,
    vector<pair<Vec3,double> >,
    int,
    int,
    int
  );
  ~GMRESSolverMaster();
};

#endif  //__GMRESSOLVERMASTER_H






