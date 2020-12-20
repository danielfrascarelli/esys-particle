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

#ifndef __SPHEREBODYINTERACTION_HPP
#define __SPHEREBODYINTERACTION_HPP

template <class T>
ASphereBodyInteraction<T>::ASphereBodyInteraction(T* p,CSphereBody* w,bool iflag)
{
  m_p=p;
  m_sphere=w;
  m_inner_flag=iflag;
  m_init=true;
}

/*!
  check if any of the particles in the interaction fits tag & mask

  \param tag the tag
  \param mask the mask
*/
template <class T>
bool ASphereBodyInteraction<T>::hasTag(int tag ,int mask) const
{
  int tag1=m_p->getTag();

  return ((tag1 & mask)==(tag & mask));
}

#endif // __SPHEREBODYINTERACTION_HPP
