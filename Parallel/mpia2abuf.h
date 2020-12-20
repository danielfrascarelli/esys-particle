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

#ifndef __MPIA2ABUFFER_H
#define __MPIA2ABUFFER_H

#include <mpi.h>

/*!
  \class CMPIA2ABuffer
  \brief class for a MPI-buffer supporting all-to-all communication
  
  \author Steffen Abe
  $Revision$
  $Date$
*/
class CMPIA2ABuffer
{
 private:
  MPI_Comm m_comm;
  int m_rank; //!< the rank in this communicator
  int m_size; //!< size of the communicator
  int m_int_increment,m_dbl_increment; //!< the "packing size" of int/double
  MPI_Status m_status;
  char* m_buffer_s; //!< send buffer
  char* m_buffer_r; //!< receive buffer
  int m_buffersize; //!< the size of the buffer per slice
  int *m_position_s; //!< the current end of the content in each slice of the send buffer
  int *m_position_r; //!< the current end of the content in each slice of the receive buffer 

 public:
  CMPIA2ABuffer(MPI_Comm,int);
  virtual ~CMPIA2ABuffer();

  virtual void clear();
  virtual void all2all();
  virtual void append(int,int);
  virtual void append(double,int);
  virtual int pop_int(int);
  virtual double pop_double(int);
};

#endif //__MPIA2ABUFFER_H
