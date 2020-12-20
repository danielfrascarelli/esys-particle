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

template< typename T>
RingBuffer<T>::RingBuffer(int s)
{
  m_buffer=vector<T>(s);
  m_idx=0;
  m_size=s;
}

template< typename T>
T& RingBuffer<T>::operator[](int i)
{
  int real_idx=(m_idx+i)%m_size;
  return m_buffer[real_idx];
}

template< typename T>
T RingBuffer<T>::operator[] (int i) const
{
  int real_idx=(m_idx+i)%m_size;
  return m_buffer[real_idx];
}

template< typename T>
void RingBuffer<T>::insert(const T& data)
{
  m_idx=(m_idx+1)%m_size;
  m_buffer[m_idx]=data;
}
