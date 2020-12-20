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

#ifndef __MPISGBUF_H
#define __MPISGBUF_H

#include <mpi.h>
#include <string>
#include "Parallel/mpibuf.h"

/*!
  \class AMPISGBufferRoot
  \brief Abstract base class for scatter/gather buffer, root component

  \author Steffen Abe
  $Revision$
  $Date$
*/
class AMPISGBufferRoot
{
protected:
  MPI_Comm m_comm; //!< the MPI communicator used for the scatter/gather operations
  int m_rank; //!< the rank in this communicator
  int m_size; //!< size of the communicator
  int m_int_increment,m_dbl_increment; //!< the "packing size" of int/double
  MPI_Status m_status;
  
public:
  AMPISGBufferRoot(MPI_Comm);
  virtual ~AMPISGBufferRoot(){};

  virtual void clear()=0;
  virtual void gather()=0;
  virtual void scatter()=0;
  virtual void append(int,int)=0;
  virtual void append(double,int)=0;
  virtual void append(const char*,int)=0;
  virtual void append(const Vec3 &,int);
  virtual int pop_int(int)=0;
  virtual double pop_double(int)=0;
  virtual void pop_doubles(int,double *,int)=0;
  virtual Vec3 pop_vector(int);
  const MPI_Status& status(){return m_status;};  
};

/*!
  \class AMPISGBufferLeaf
  \brief Abstract base class for scatter/gather buffer, leaf component

  \author Steffen Abe
  $Revision$
  $Date$
*/
class AMPISGBufferLeaf : public AMPIBuffer
{
 protected:
  int m_root; //!< rank of the root process
  int m_int_increment,m_dbl_increment; //!< the "packing size" of int/double

public:
  AMPISGBufferLeaf(MPI_Comm,int);
  virtual ~AMPISGBufferLeaf(){};

  virtual void clear()=0;
  virtual void send()=0;
  virtual void receive()=0;
  virtual void append(int)=0;
  virtual void append(double)=0;
  virtual int pop_int()=0;
  virtual double pop_double()=0;
  virtual void pop_doubles(double *,int)=0;
  virtual std::string pop_string()=0;
  const MPI_Status& status(){return m_status;};  
};

/*!
  \class CMPISGBufferRoot
  \brief Buffer for MPI scatter/gather, root component

  \author Steffen Abe
  $Revision$
  $Date$
*/
class CMPISGBufferRoot : public AMPISGBufferRoot
{
private:
  char* m_buffer;
  char* m_dummy_buffer; //!<dummy buffer sent by root to itself
  int m_buffersize; //!< the size of the buffer per slice
  int *m_position; //!< the current end of the content in each slice

public:
  CMPISGBufferRoot(MPI_Comm,int);
  virtual ~CMPISGBufferRoot();

  virtual void clear();
  virtual void gather();
  virtual void scatter();
  virtual void append(int,int);
  virtual void append(double,int);
  virtual void append(const char*,int);
  virtual int pop_int(int);
  virtual double pop_double(int);
  virtual void pop_doubles(int,double *,int);
};

/*!
  \class CMPISGBufferLeaf
  \brief Buffer for MPI scatter/gather, leaf component

  \author Steffen Abe
  $Revision$
  $Date$
*/
class CMPISGBufferLeaf : public AMPISGBufferLeaf
{
private:
  char* m_buffer; 
  int m_buffersize; //!< the size of the buffer
  int m_position; //!< the current end of the content

public:
  CMPISGBufferLeaf(MPI_Comm,int,int);
  virtual ~CMPISGBufferLeaf();

  virtual void clear();
  virtual void send();
  virtual void receive();
  virtual void append(int);
  virtual void append(double);
  virtual void append(const char*);
  virtual int pop_int();
  virtual double pop_double();
  virtual void pop_doubles(double *,int);
  virtual std::string pop_string();
};

#endif // __MPISGBUF_H
