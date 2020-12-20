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

#ifndef __HANDLE_H
#define __HANDLE_H

//--- project includes ---
#include "handle_exception.h"

/*!
  \class T_Handle
  \brief Template class for a handle/ref. counted pointer

  Extended version of the example in Stroustrup, p. 783
*/
template <typename T>
class T_Handle
{
 private:
  T* m_rep;
  int* m_count;
  
 public:
  T_Handle(T*);
  T_Handle(const T_Handle&);
  ~T_Handle();

  T_Handle& operator=(const T_Handle&);
  inline T* operator->() {return m_rep;};
  inline T& operator*(){return *m_rep;};
  //void destroy();
};

#include "handle.hpp"

#endif // __HANDLE_H
