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

#ifndef __PACKED_MESSAGE_INTERFACE_H
#define __PACKED_MESSAGE_INTERFACE_H

//--- project includes ---
#include "Foundation/vec3.h" // for append(Vec3), pop_vec3()
#include "Foundation/Matrix3.h"

//--- STL includes ---
#include <string>
using std::string;

/*!
  \class TML_PackedMessageInterface
  \brief Abstract base/interface class for packed messages to be used in
  TML_Pack
*/
class TML_PackedMessageInterface
{
 public:
  virtual void begin_pack()=0;
  virtual void begin_unpack()=0;
  virtual void append(int)=0; 
  virtual void append(double)=0; 
  virtual void append(const string&)=0;
  virtual void append(const Vec3&)=0;
//  virtual void append(const Matrix3&)=0;
  virtual void append(bool)=0;

  virtual int pop_int()=0;
  virtual double pop_double()=0;
  virtual void pop_doubles(double*,int)=0;
  virtual string pop_string()=0;
  virtual Vec3 pop_vec3()=0;
 // virtual Matrix3 pop_matrix3()=0;
  virtual bool pop_bool()=0;

  template<typename T> void pack(const T&);
  template<typename T> void unpack(T&);
};

#endif //__PACKED_MESSAGE_INTERFACE_H
