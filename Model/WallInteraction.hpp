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

#ifndef __WALLINTERACTION_HPP
#define __WALLINTERACTION_HPP

template <class T>
AWallInteraction<T>::AWallInteraction(T* p,CWall* w,bool iflag)
{
  m_p=p;
  m_wall=w;
  m_inner_flag=iflag;
  m_init=true;
}

/*!
  check if any of the particles in the interaction fits tag & mask

  \param tag the tag
  \param mask the mask
*/
template <class T>
bool AWallInteraction<T>::hasTag(int tag ,int mask) const
{
  int tag1=m_p->getTag();

  return ((tag1 & mask)==(tag & mask));
}

#endif  // for the error of redefination ! 
