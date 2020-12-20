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

#include "mpicmdbuf.h"
#include "console.h"

/*!
  Constructor. Sets the MPI communicator to be used for broadcast operations and the rank of root process;

  \param comm the MPI communicator
  \param root the rank of the root process

*/
CMPILCmdBuffer::CMPILCmdBuffer(MPI_Comm comm,int root)
{
  int rank;

  m_comm=comm;
  m_root=root;
  MPI_Comm_rank(m_comm, &rank);
  m_isroot=(rank==root);
}

/*!
  Broadcast a command to all members of the communicator. If the calling process is not the root process, prints an error message and does nothing.

  \param cmd the command
*/
void CMPILCmdBuffer::broadcast(int cmd)
{
  if(m_isroot){
    MPI_Bcast(&cmd,1,MPI_INT,m_root,m_comm);
  } else {
    console.Error() << "Error in CMPILCmdBuffer::broadcast : trying to broadcast from non-root process !\n";
  }
}

/*!
  receive broadcast and return the received command
*/
int CMPILCmdBuffer::receive()
{
  int cmd;
  MPI_Bcast(&cmd,1,MPI_INT,m_root,m_comm);
  return cmd;
}
