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

#ifndef __MPIVBUF_H
#define __MPIVBUF_H

#include <mpi.h>
#include <string>
#include "Parallel/mpibuf.h"

/*!
  \class CVarMPIBuffer
  \brief MPI send/recv buffer with automagically adjusted size
  
  CVarMPIBuffer implements a send/receive buffer with variable size. The buffer grows automatically if neccesary with each append operation and when the a message is received. It is never automatically shrunk. Both append and send/recv operations slower than a constant size buffer (CMPIBuffer)  

  \author Steffen Abe
  $Revision$
  $Date$

  \todo implement checks for locking
  \todo use exeption handling for error checking
*/

class CVarMPIBuffer : public AMPIBufferPP
{
 private:
  char* m_buffer;
  int m_buffersize; //!< the size of the buffer
  int m_position; //!< the current end of the content
  int m_int_increment,m_dbl_increment; //!< the "packing size" of int/double
  bool m_lock; 
  MPI_Request m_req[2];
  MPI_Status m_stat[2];//!< for the ISendTo/Wait stuff

 protected:
  void grow();
  void growTo(int);

 public:
  CVarMPIBuffer(MPI_Comm,int size=16);
  virtual ~CVarMPIBuffer();

  virtual void clear(){m_position=0;};
  virtual void sendTo(int,int);
  virtual void NBsendTo(int,int);
  virtual void initSendTo(int,int);
  virtual void wait();
  virtual void receiveFrom(int src=MPI_ANY_SOURCE,int tag=MPI_ANY_TAG);
  virtual void append(int);
  virtual void append(double);
  virtual void append(const char*);
  virtual void append(const Vec3 &V) { AMPIBufferPP::append(V); } ;
  virtual int pop_int();
  virtual double pop_double();
  virtual void pop_doubles(double*,int);
  virtual std::string pop_string();
  virtual void broadcast(int);
  virtual void receiveBroadcast(int);
};
#endif //__MPIVBUF_H
