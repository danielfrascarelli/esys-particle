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


#include "tml/message/multi_message_slab.h"
#include "tml/message/packed_multi_message.h"

/*!
  construct TML_PackedMultiMessageSlab

  \param msg the multimessage to which the slab belongs
  \param idx the index 
*/
TML_PackedMultiMessageSlab::TML_PackedMultiMessageSlab(TML_PackedMultiMessage* msg,int idx)
{
  m_msg=msg;
  m_idx=idx;
}

/*!
  reset the packing pointer
*/
void TML_PackedMultiMessageSlab::begin_pack()
{
  m_msg->begin_pack(m_idx);
}

/*!
  reset the unpacking pointer
*/
void TML_PackedMultiMessageSlab::begin_unpack()
{
  m_msg->begin_unpack(m_idx);
}
  
/*!
  pack integer into the slab

  \param i the interger
*/
void TML_PackedMultiMessageSlab::append(int i)
{
  m_msg->append(i,m_idx);
}
 
/*!
  pack a double into the slab

  \param d the double
*/
void TML_PackedMultiMessageSlab::append(double d)
{
  m_msg->append(d,m_idx);
}

/*!
  pack a STL string into the slab

  \param s the string
*/
void TML_PackedMultiMessageSlab::append(const string& s)
{
  m_msg->append(s,m_idx);
}

/*!
  Append a Vec3 to the message buffer. Calls append(double) per element
*/
void TML_PackedMultiMessageSlab::append(const Vec3& v)
{
  m_msg->append(v[0],m_idx);
  m_msg->append(v[1],m_idx);
  m_msg->append(v[2],m_idx);
}

/*!
  pack a booleam value into the slab

  \param b the boolean value
*/
void TML_PackedMultiMessageSlab::append(bool b)
{
  m_msg->append(b,m_idx);
}

/*!
  pop an int from the slab
*/
int TML_PackedMultiMessageSlab::pop_int()
{
  return m_msg->pop_int(m_idx);
}

/*!
  pop a double from the slab
*/
double TML_PackedMultiMessageSlab::pop_double()
{
  return m_msg->pop_double(m_idx);
}

/*!
  pop an array of doubles from the slab
  \warning not implemented
*/
void TML_PackedMultiMessageSlab::pop_doubles(double*,int)
{}

/*!
  pop a STL string  from the slab
  \warning not implemented
*/
string TML_PackedMultiMessageSlab::pop_string()
{
  string s;
  return s;
}

/*!
  Pop a Vec3 of the buffer. Calls pop_double per element
*/
Vec3 TML_PackedMultiMessageSlab::pop_vec3()
{
  Vec3 res;

  res[0]=m_msg->pop_double(m_idx);
  res[1]=m_msg->pop_double(m_idx);
  res[2]=m_msg->pop_double(m_idx);

  return res;
}

/*!
  pop a boolean value from the slab
*/
bool TML_PackedMultiMessageSlab::pop_bool()
{
  return m_msg->pop_bool(m_idx);
}
