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

#include "BasicParticle.h"

CBasicParticle::CBasicParticle()
{
  m_id=-1;
  m_pos=Vec3(0.0,0.0,0.0);
  m_rad=0.0;
}

CBasicParticle::CBasicParticle(int id,const Vec3& pos,double rad)
{
  m_id=id;
  m_pos=pos;
  m_rad=rad;
}
