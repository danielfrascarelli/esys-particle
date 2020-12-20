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

#ifndef __MPISGVBUF_H
#define __MPISGVBUF_H

#include "Parallel/mpisgbuf.h"
#include <string>

/*!
  \class CMPIVarSGBufferRoot
  \brief class for variable size scatter/gather buffer, root component
  
  \author Steffen Abe
  $Revision$
  $Date$
*/
class CMPIVarSGBufferRoot: public AMPISGBufferRoot
{ 
private:
  char* m_vbuffer;
  char* m_dummy_vbuffer; //!<dummy buffer sent by root to itself
  int m_vbuffersize; //!< the size of the buffer per slice
  int *m_position; //!< the current end of the content in each slice
  int *m_rpos;     //!< the number of bytes in the slice (i.e. m_position-m_displ)

  int *m_recvcount;//!< the buffer for the transfer of the size of the vbuffer
  int *m_displ; //<! the diplacements of the slices in the buffer
  int m_ndummy;
 

 protected:
  void grow();
  void growTo(int);

 public:
  CMPIVarSGBufferRoot(MPI_Comm,int isize=16);
  virtual ~CMPIVarSGBufferRoot();

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
  \class CMPIVarSGBufferLeaf
  \brief class for variable size scatter/gather buffer, leaf component
  
  \author Steffen Abe
  $Revision$
  $Date$
*/
class CMPIVarSGBufferLeaf: public AMPISGBufferLeaf
{
 private:
  char* m_vbuffer; 
  int m_vbuffersize; //!< the size of the buffer
  int m_position; //!< the current end of the content
  int m_data_size;

 protected:
  void grow();
  void growTo(int);

 public:
  CMPIVarSGBufferLeaf(MPI_Comm,int,int isize=16);
  virtual ~CMPIVarSGBufferLeaf();

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



#endif // __MPISGVBUF_H
