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

#ifndef __PACKED_MESSAGE_H
#define __PACKED_MESSAGE_H

//--- MPI ---
#include <mpi.h>

//--- project includes ---
#include "tml/message/packed_message_interface.h"
#include "Foundation/vec3.h" // for append(Vec3), pop_vec3()
#include "Foundation/Matrix3.h"

//--- STL includes ---
#include <string>
using std::string;


/*!
  \class TML_Packed_Message
  \brief Message buffer for sending and receiving packed data. Data types are not checked.  
  The implementation has been derived from the old CVarMPIBuffer class 
  
*/
class TML_Packed_Message : public TML_PackedMessageInterface
{
 protected:
  MPI_Comm m_comm;  
  char* m_buffer;
  int m_buffersize; //!< the size of the buffer
  int m_pack_pos; //!< the current end of the content
  int m_unpack_pos; //!< the current pos for unpacking
  int m_int_increment,m_dbl_increment; //!< the "packing size" of int/double

  void grow();
  void growTo(int);

 public:
  TML_Packed_Message(MPI_Comm,unsigned int size=64);
  virtual ~TML_Packed_Message();

  char* buffer(){return m_buffer;}; // make protected ?
  int size(){return m_pack_pos;};

  virtual void begin_pack(){m_pack_pos=0;};
  virtual void begin_unpack(){m_unpack_pos=0;};
  virtual void append(int); 
  virtual void append(double); 
  virtual void append(const string&);
  virtual void append(const Vec3&);
  virtual void append(const Matrix3&);
  virtual void append(bool);

  virtual int pop_int();
  virtual double pop_double();
  virtual void pop_doubles(double*,int);
  virtual string pop_string();
  virtual Vec3 pop_vec3();
  virtual Matrix3 pop_matrix3();
  virtual bool pop_bool();
};
#endif //__PACKED_MESSAGE_H
