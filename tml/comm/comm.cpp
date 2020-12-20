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


#include "tml/comm/comm.h"
#include "Foundation/console.h"

//--- I/O ---

/*!
  Default constructor for TML_comm. Sets the MPI communicator
  to MPI_COMM_NULL;
*/
TML_Comm::TML_Comm()  
{
  m_comm=MPI_COMM_NULL;
}

/*!
  construct TML_Comm from MPI communicator

  \param comm the MPI communicator
*/
TML_Comm::TML_Comm(MPI_Comm comm)
{
  m_comm=comm;
}

/*!
  set the undelying MPI communicator

  \param comm the MPI communicator
*/
void TML_Comm::setComm(MPI_Comm comm)
{
  m_comm=comm;
}

TML_Comm& TML_Comm::operator=(const TML_Comm& rhs)
{
  //  MPI_Comm_dup(rhs.m_comm,&m_comm);
  m_comm=rhs.m_comm;

  return *this;
}

int TML_Comm::rank()  const
{
  int rank;
  
  if(m_comm!=MPI_COMM_NULL){
    MPI_Comm_rank(m_comm,&rank);
  } else {
    rank=MPI_UNDEFINED;
  }

  return rank;
}

int TML_Comm::size()
{
  int size;

  if(m_comm!=MPI_COMM_NULL){
    MPI_Comm_size(m_comm,&size);
  } else {
    size=MPI_UNDEFINED;
  }

  return size;
}

/*!
  Construct a new communicator containing the processes
  which are given as input

  \param ids the ranks of the processes which form the new communicator in the current communicator
  \todo error handling
*/
TML_Comm TML_Comm::include(const vector<int>& ids)
{
  TML_Comm newcomm;
  MPI_Group grp,ngrp;

  // extract group
  MPI_Comm_group(m_comm,&grp);
  // vector->array
  int nids=ids.size();
  int *ranks=new int[nids];
  for(int i=0;i<nids;i++){
    ranks[i]=ids[i];
  }
  // make new group
  int err=MPI_Group_incl(grp,nids,ranks,&ngrp);
  if(err!=MPI_SUCCESS){
    console.Error() << "Error in TML_Comm::include group construction, rank " 
	 << rank() << " error " << err << "\n";
  }
  int gsize,grnk;
  MPI_Group_size(ngrp,&gsize);
  MPI_Group_rank(ngrp,&grnk);
  delete ranks;
  // construct new MPI communicator
  err=MPI_Comm_create(m_comm,ngrp,&(newcomm.m_comm));

  if(err!=MPI_SUCCESS){
    console.Error() << "Error in TML_Comm::include communicator construction, rank " 
	 << rank() << " error " << err << "\n";
  }

  return newcomm;
}

/*!
  Construct a new communicator containing the processes from the current communicator 
  except to ones which are given as input

  \param ids the ranks of the processes which are excluded from  the new communicator 
  \todo error handling
*/
TML_Comm TML_Comm::exclude(const vector<int>& ids)
{
  TML_Comm newcomm;
  MPI_Group grp,ngrp;

  // extract group
  MPI_Comm_group(m_comm,&grp);
  // vector->array
  int nids=ids.size();
  int *ranks=new int[nids];
  for(int i=0;i<nids;i++){
    ranks[i]=ids[i];
  }
  // make new group
  MPI_Group_excl(grp,nids,ranks,&ngrp);
  delete ranks;
  // construct new MPI communicator
  MPI_Comm_create(m_comm,ngrp,&(newcomm.m_comm));

  return newcomm;
}

/*!
  Wait on a barrier. Wrapper for MPI_Barrier.
*/
void TML_Comm::barrier()
{
  MPI_Barrier(m_comm);
}

/*!
  Wait on a barrier with debug message
*/
void TML_Comm::barrier(const string& msg)
{
  double m_time=MPI_Wtime();
  int id=rank();
  if(id==0){
    console.Debug() << "Master waiting on Barrier ( " << msg << " )\n";
  } else {
    console.Debug() << "Worker " << id << " waiting on Barrier ( " << msg << " )\n";
  }

  MPI_Barrier(m_comm);
  double p_time=MPI_Wtime();
  if(id==0){
    console.Debug() << "Master past Barrier ( " << msg << " ) after " << p_time-m_time << " sec \n";
  } else {
    console.Debug() << "Worker " << id << " past Barrier ( " << msg << " ) after " << p_time-m_time << " sec \n";
  }
}
