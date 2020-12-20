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

#ifndef __MULTI_MESSAGE_SLAB_H
#define __MULTI_MESSAGE_SLAB_H

//--- MPI ---
#include <mpi.h>


//--- project includes---
#include "tml/message/packed_message_interface.h"
#include "Foundation/vec3.h" // for append(Vec3), pop_vec3()

//--- forward declarations ---

class TML_PackedMultiMessage;

/*!
  \class TML_PackedMultiMessageSlab
  \brief Handle class to access multimessages via a packed message interface
*/
class TML_PackedMultiMessageSlab : public TML_PackedMessageInterface
{
 private:
  TML_PackedMultiMessage *m_msg;
  int m_idx;

 public:
  TML_PackedMultiMessageSlab(TML_PackedMultiMessage*,int);

  virtual void begin_pack();
  virtual void begin_unpack();
  virtual void append(int); 
  virtual void append(double); 
  virtual void append(const string&);
  virtual void append(const Vec3&);
  virtual void append(bool);

  virtual int pop_int();
  virtual double pop_double();
  virtual void pop_doubles(double*,int);
  virtual string pop_string();
  virtual Vec3 pop_vec3();
  virtual bool pop_bool();
};
#endif //__MULTI_MESSAGE_SLAB_H
