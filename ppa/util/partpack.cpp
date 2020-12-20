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

#include "packed_message_interface.h"
#include "BasicParticle.h"
#include "vec3.h"

template<>
void TML_PackedMessageInterface::pack<CBasicParticle>(const CBasicParticle& p)
{
  Vec3 pos=p.getPos();
  append(pos.X());
  append(pos.Y());
  append(pos.Z());
  append(p.getRad());
  append(p.getID());
}

template<>
void TML_PackedMessageInterface::unpack<CBasicParticle>(CBasicParticle& p)
{
  double x=pop_double();
  double y=pop_double();
  double z=pop_double();
  double r=pop_double();
  int id=pop_int();
  p=CBasicParticle(id,Vec3(x,y,z),r);
}
