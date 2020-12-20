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

#ifndef __MPIBUF_H
#define __MPIBUF_H

//--- MPI includes ---
#include <mpi.h>

//--- Project includes ---
#include "Foundation/vec3.h"

// --- STL includes ---
#include <string>

/*!
  \class AMPIBuffer
  \brief Abstract base class for MPI send/recv buffer

  \author Steffen Abe
  $Revision$
  $Date$
*/
class AMPIBuffer
{
 protected:
  MPI_Comm   m_comm; //!< the MPI Communicator used for the send/recv operations
  MPI_Status m_status;

public:
  AMPIBuffer(MPI_Comm comm){m_comm=comm;};
  virtual ~AMPIBuffer(){};

  virtual void clear()=0;
  virtual void append(int)=0;
  virtual void append(double)=0;
  virtual void append(const char*)=0;
  virtual void append(const Vec3 &) ;
  virtual int pop_int()=0;
  virtual double pop_double()=0;
  virtual void pop_doubles(double*,int)=0;
  virtual std::string pop_string()=0;
  virtual Vec3 pop_vector() ;
  const MPI_Status& status(){return m_status;};   
};

/*!
  \class AMPIBufferPP
  \brief Abstarct base class for Point-to-Point communication buffers

  Adds sendTo and receiveFrom to the base class
  \author Steffen Abe
  $Revision$
  $Date$
*/
class AMPIBufferPP : public AMPIBuffer
{
public:
  AMPIBufferPP(MPI_Comm comm);
  virtual ~AMPIBufferPP(){};
  virtual void sendTo(int,int)=0;
  virtual void receiveFrom(int src=MPI_ANY_SOURCE,int tag=MPI_ANY_TAG)=0;
};

/*!
  \class CMPIBuffer
  \brief Constant size MPI send/recv buffer

  CMPIBuffer implements a send/receive buffer. MPI_pack/MPI_unpack is used the transfer arbitrary data. Type information is not transported. i.e. the user has to know the type of the content of a received message.

  \author Steffen Abe
  $Revision$
  $Date$
*/

class CMPIBuffer : public AMPIBufferPP
{
 private:
 
  char* m_buffer;
  int m_buffersize; //!< the size of the buffer
  int m_position; //!< the current end of the content
  int m_int_increment,m_dbl_increment; //!< the "packing size" of int/double
  

 public:
  CMPIBuffer(MPI_Comm,int);
  virtual ~CMPIBuffer();

  virtual void clear(){m_position=0;};
  virtual void sendTo(int,int);
  virtual void receiveFrom(int src=MPI_ANY_SOURCE,int tag=MPI_ANY_TAG);
  virtual void append(int);
  virtual void append(double);
  virtual void append(const char*);
  bool append_checked(int);
  bool append_checked(double);
  virtual int pop_int();
  virtual double pop_double();
  virtual void pop_doubles(double*,int);
  virtual std::string pop_string();
};

#endif //__MPIBUF_H
