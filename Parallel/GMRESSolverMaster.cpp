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


#include "GMRESSolverMaster.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "Parallel/mpicmdbuf.h"
#include "Parallel/mpisgvbuf.h"
#include "Parallel/sublattice_cmd.h"
#include "Parallel/mpibarrier.h"
#include "Parallel/BroadCast_cmd.h"
#include "Parallel/mpi_tag_defs.h"
#include "tml/comm/comm.h"

using namespace std;


GMRESSolverMaster::GMRESSolverMaster(
  MPI_Comm* comm,
  TML_Comm* tml_comm,
  vector<pair<Vec3,double> > co_W,
  vector<pair<Vec3,double> > co_E,
  vector<pair<Vec3,double> > co_N,
  vector<pair<Vec3,double> > co_S,
  vector<pair<Vec3,double> > co_D,
  vector<pair<Vec3,double> > co_U,
  vector<pair<Vec3,double> > co_C,
  vector<pair<Vec3,double> > co_B,
  int nx,
  int ny,
  int nz
)
{
  m_comm=comm;
  m_tml_comm=tml_comm;
  MPI_Comm_rank(*m_comm, &m_rank);

  Nx=nx;Ny=ny;Nz=nz;
  Number=Nx*Ny*Nz;
  W.resize(Number); E.resize(Number); N.resize(Number); S.resize(Number); D.resize(Number); U.resize(Number); C.resize(Number); B.resize(Number);
  irow.resize(Number*7); icol.resize(Number*7); var.resize(Number*7);
  index.resize(Number);
  for(int n=0;n<Number;n++){
    W[n]=(&co_W[n])->second; E[n]=(&co_E[n])->second;
    N[n]=(&co_N[n])->second; S[n]=(&co_S[n])->second;
    D[n]=(&co_D[n])->second; U[n]=(&co_U[n])->second;
    C[n]=(&co_C[n])->second; B[n]=(&co_B[n])->second;
    index[n]=(&co_W[n])->first;
  }
}


GMRESSolverMaster::~GMRESSolverMaster()
{
  irow.clear();icol.clear();
  W.clear();E.clear();N.clear();S.clear();D.clear();U.clear();C.clear();B.clear();
  var.clear();index.clear();
}


vector<pair<Vec3,double> > GMRESSolverMaster::MatrixSolving()
{
  CMPILCmdBuffer cmd_buffer(*m_comm,m_rank);
  CMPIBarrier barrier(*m_comm);

  GlobalMatrixSetting();

  cmd_buffer.broadcast(CMD_SOLVE);

  m_tml_comm->broadcast_cont(B);
  m_tml_comm->broadcast_cont(irow);
  m_tml_comm->broadcast_cont(icol);
  m_tml_comm->broadcast_cont(var);
  m_tml_comm->broadcast(Number);
  m_tml_comm->broadcast(nnz);

  multimap<int,double> temp;
  m_tml_comm->gather(temp);

  barrier.wait("solve matrix");

  vector<double> x;
  for(multimap<int,double>::iterator iter=temp.begin();
      iter!=temp.end();
      iter++){
    x.push_back(iter->second);
  };

  //putting x into results vector
  vector<pair<Vec3,double> > results;
  for(int n=0;n<Number;n++){
    results.push_back(make_pair(index[n],x[n]));
  }

  x.clear();
  temp.clear();
  return results;
}


void GMRESSolverMaster::GlobalMatrixSetting()
{
  nnz=0;
  for(int n=1,m1=0,m2=0;n<=Number;n++){
    if (n>Nx*Ny){irow[nnz]=n;icol[nnz]=n-Nx*Ny;var[nnz]=D[n-1];nnz++;};//down3
    if (n>Nx*(1+Ny*m1) && n<=Nx*Ny*(m1+1)) {
      irow[nnz]=n;icol[nnz]=n-Nx;var[nnz]=N[n-1];nnz++; //down2
      if(n==Nx*Ny*(m1+1)){m1++;};
    };
    if(n%Nx!=1){irow[nnz]=n;icol[nnz]=n-1;var[nnz]=W[n-1];nnz++;};//down1
    irow[nnz]=n;icol[nnz]=n;var[nnz]=C[n-1];nnz++; //centre
    if(n%Nx!=0){irow[nnz]=n;icol[nnz]=n+1;var[nnz]=E[n-1];nnz++;};//up1
    if (n>m2*Nx*Ny && n<=(m2+1)*Ny*Nx-Nx) {
      irow[nnz]=n;icol[nnz]=n+Nx;var[nnz]=S[n-1];nnz++;//up2
      if(n==(m2+1)*Ny*Nx-Nx){m2++;};
    };
    if (n<=Nx*Ny*(Nz-1)){irow[nnz]=n;icol[nnz]=n+Nx*Ny;var[nnz]=U[n-1];nnz++;};//up3
  };
}



