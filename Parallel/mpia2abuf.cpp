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

#include "mpia2abuf.h"

CMPIA2ABuffer::CMPIA2ABuffer(MPI_Comm comm,int buffersize)
{
  m_comm=comm;
  MPI_Comm_rank(comm,&m_rank);
  MPI_Comm_size(comm,&m_size);
  MPI_Pack_size(1,MPI_INT,m_comm,&m_int_increment);
  MPI_Pack_size(1,MPI_DOUBLE,m_comm,&m_dbl_increment);

  m_buffersize=buffersize;
  m_buffer_s=new char[m_buffersize*m_size];
  m_buffer_r=new char[m_buffersize*m_size];
  m_position_s=new int[m_size];
  m_position_r=new int[m_size];
  for(int i=0;i<m_size;i++){
    m_position_s[i]=i*m_buffersize;
    m_position_r[i]=i*m_buffersize;
  }
}

CMPIA2ABuffer::~CMPIA2ABuffer()
{
  delete m_buffer_s;
  delete m_buffer_r;
  delete m_position_s;
  delete m_position_r;
}

void CMPIA2ABuffer::clear()
{
  for(int i=0;i<m_size;i++){
    m_position_s[i]=i*m_buffersize;
    m_position_r[i]=i*m_buffersize;
  } 
}

void CMPIA2ABuffer::all2all()
{
  MPI_Alltoall(m_buffer_s,m_buffersize,MPI_PACKED,m_buffer_r,m_buffersize,MPI_PACKED,m_comm);
}

/*!
  Append an integer to a given slice of the buffer.

  \param i the integer
  \param nslice the nr. of the slice

  \warning No check for overflow
*/
void CMPIA2ABuffer::append(int i,int nslice)
{
  MPI_Pack(&i,1,MPI_INT,m_buffer_s,m_buffersize*m_size,&(m_position_s[nslice]),m_comm);  
}

/*!
  Append an double to a given slice of the buffer.

  \param d the double
  \param nslice the nr. of the slice

  \warning No check for overflow
*/
void CMPIA2ABuffer::append(double d,int nslice)
{
  MPI_Pack(&d,1,MPI_DOUBLE,m_buffer_s,m_buffersize*m_size,&(m_position_s[nslice]),m_comm);
}

/*!
  Pops an integer from a given slice of the the buffer

  \param nslice the nr. of the slice
  \return the int.

  \warning No check for underflow
*/
int CMPIA2ABuffer::pop_int(int nslice)
{
   int res;
  MPI_Unpack(m_buffer_r,m_buffersize*m_size,&(m_position_r[nslice]),&res,1,MPI_INT,m_comm);
  
  return res;  
}

/*!
  Pops an double from a given slice of the the buffer. 

  \param nslice the nr. of the slice
  \return the double.

  \warning No check for underflow
*/
double CMPIA2ABuffer::pop_double(int nslice)
{
  double res;
  MPI_Unpack(m_buffer_r,m_buffersize,&(m_position_r[nslice]),&res,1,MPI_DOUBLE,m_comm);

  return res;
}
