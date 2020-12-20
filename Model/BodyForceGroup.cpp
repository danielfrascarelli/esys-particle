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


#include "Model/BodyForceGroup.h"

namespace esys
{
  namespace lsm
  {
    BodyForceIGP::BodyForceIGP() : AIGParam(), m_acceleration()
    {
    }

    BodyForceIGP::BodyForceIGP(const std::string &name, const Vec3 &acceleration)
      : AIGParam(name),
        m_acceleration(acceleration)
    {
    }

    BodyForceIGP::~BodyForceIGP()
    {
    }

    const Vec3 &BodyForceIGP::getAcceleration() const
    {
      return m_acceleration;
    }

    const std::string &BodyForceIGP::getName() const
    {
      return Name();
    }

    void BodyForceIGP::packInto(CVarMPIBuffer *pBuffer) const
    {
      pBuffer->append(getName().c_str());
      pBuffer->append(m_acceleration.X());
      pBuffer->append(m_acceleration.Y());
      pBuffer->append(m_acceleration.Z());
    }

    BodyForceIGP BodyForceIGP::extract(CVarMPIBuffer *pBuffer)
    {
      const std::string name = pBuffer->pop_string();

      const double x = pBuffer->pop_double();
      const double y = pBuffer->pop_double();
      const double z = pBuffer->pop_double();

      return BodyForceIGP(name, Vec3(x, y, z));
    }

    BuoyancyIGP::BuoyancyIGP (const std::string &name, const Vec3 &acceleration, const double &fluidDensity, const double &fluidHeight) : AIGParam(name), m_acceleration(acceleration), m_fluidDensity (fluidDensity), m_fluidHeight (fluidHeight)
    {
    }

    const Vec3 &BuoyancyIGP::getAcceleration() const
    {
      return m_acceleration;
    }

    const double &BuoyancyIGP::getFluidDensity() const
    {
      return m_fluidDensity;
    }

    const double &BuoyancyIGP::getFluidHeight() const
    {
      return m_fluidHeight;
    }

    void BuoyancyIGP::packInto(CVarMPIBuffer *pBuffer) const
    {
      pBuffer->append(getName().c_str());
      pBuffer->append(m_acceleration.X());
      pBuffer->append(m_acceleration.Y());
      pBuffer->append(m_acceleration.Z());
      pBuffer->append(m_fluidDensity);
      pBuffer->append(m_fluidHeight);
    }

    BuoyancyIGP BuoyancyIGP::extract(CVarMPIBuffer *pBuffer)
    {
      const std::string name = pBuffer->pop_string();

      const double x = pBuffer->pop_double();
      const double y = pBuffer->pop_double();
      const double z = pBuffer->pop_double();
      const double fluidDensity = pBuffer->pop_double();
      const double fluidHeight = pBuffer->pop_double();

      return BuoyancyIGP(name, Vec3(x, y, z), fluidDensity, fluidHeight);
    }
  }
}
