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

#ifndef __RINGBUFFER_H
#define __RINGBUFFER_H

// --- STL includes ---
#include <vector>
using std::vector;

/*!
  class for a ringbuffer
*/
template <typename T>
class RingBuffer
{
 private:
  vector<T> m_buffer;
  int m_idx;
  int m_size;

 public:
  RingBuffer(int);

  T& operator[](int);
  T operator[] (int) const;
  void insert(const T&);
  inline int size() const {return m_size;};
};

#include "Foundation/RingBuffer.hpp"

#endif // __RINGBUFFER_H
