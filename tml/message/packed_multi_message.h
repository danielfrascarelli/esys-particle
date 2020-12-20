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

#ifndef __PACKED_MULTI_MESSAGE_H
#define __PACKED_MULTI_MESSAGE_H

//--- MPI ---
#include <mpi.h>

//--- project includes ---
#include "tml/message/multi_message_slab.h"

//--- STL includes ---
#include <string>
using std::string;


/*!
  \class TML_PackedMultiMessage
  \brief Message buffer for sending and receiving packed data to mutltiple receivers. 
  Data types are not checked. The implementation has been derived from the old 
  CMPIVarSGBufferRoot class
  
*/
class TML_PackedMultiMessage
{
 private:
  MPI_Comm m_comm;  
  char* m_vbuffer;
  int m_vbuffersize; //!< the size of the buffer per slice
  int *m_position; //!< the current end of the content in each slice
  int *m_rpos;     //!< the number of bytes in the slice (i.e. m_position-m_displ)

  int *m_recvcount;//!< the buffer for the transfer of the size of the vbuffer
  int *m_displ; //<! the diplacements of the slices in the buffer
  int m_size;
  int m_int_increment,m_dbl_increment; //!< the "packing size" of int/double

 protected:
  void grow();
  void growTo(int);

 public:
  TML_PackedMultiMessage(MPI_Comm,int isize=64);
  virtual ~TML_PackedMultiMessage();

  TML_PackedMultiMessageSlab operator[](int);

  char* buffer(){return m_vbuffer;}; // make protected or make TML_Comm friend?
  int* offsets(){return m_displ;};
  int* sizes(){return m_rpos;};

  void clear();
  void begin_pack(int);
  void begin_unpack(int);

  void append(int,int);
  void append(double,int);
  void append(const string&,int);
  void append(bool,int);

  int pop_int(int);
  double pop_double(int);
  //string pop_string(int);
  string pop_string();
  bool pop_bool(int);
};

#endif // __PACKED_MULTI_MESSAGE_H
