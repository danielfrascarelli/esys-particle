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

#include "Model/BasicParticle.h"
#include "Geometry/SimpleParticleData.h"

const CBasicParticle CBasicParticle::INVALID = CBasicParticle(Vec3::ZERO, 0.0, -1, -1);

CBasicParticle::CBasicParticle()
  : m_pos(0.0,0.0,0.0),
    m_rad(0.0),
    m_global_id(-1),
    m_tag(-1)
{
}

CBasicParticle::CBasicParticle(const esys::lsm::SimpleParticleData &data)
  : m_pos(data.getPosition()),
    m_rad(data.getRadius()),
    m_global_id(data.getId()),
    m_tag(data.getTag())
{
}

CBasicParticle::CBasicParticle(const Vec3& pos, double rad, int id, int tag)
  : m_pos(pos),
    m_rad(rad),
    m_global_id(id),
    m_tag(-1)
{
}

ostream& operator<<(ostream& ost,const CBasicParticle& BP)
{
  ost << "Particle- id " << BP.getID() << " pos: " << BP.getPos() << " rad: " << BP.getRad() << " tag : " << BP.getTag() << endl;  
  return ost;
}
