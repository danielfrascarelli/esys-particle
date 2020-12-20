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

#ifndef ESYS_LSMCONTACT_H
#define ESYS_LSMCONTACT_H

#include "Foundation/vec3.h"

#include <math.h>
#include <iostream>

namespace esys
{
  namespace lsm
  {
    class ParticleData
    {
    public:
      ParticleData() : m_pos(), m_rad()
      {
      }

      ParticleData(const Vec3 &pos, double rad) : m_pos(pos), m_rad(rad)
      {
      }

      bool operator==(const ParticleData &pd) const
      {
        return ((getPos() == pd.getPos()) && (getRad() == pd.getRad()));
      }

      const Vec3 &getPos() const
      {
        return m_pos;
      }
      
      const double &getRad() const
      {
        return m_rad;
      }
      
      static bool is3d()
      {
        return s_is3d;
      }

      static bool is2d()
      {
        return !s_is3d;
      }

      static bool is3d(bool is3d)
      {
        s_is3d = is3d;
        return s_is3d;
      }
      
      static const double FOUR_THIRDS_PI;
      double getVolume() const
      {
        return 
          (
            is3d()
            ?
              FOUR_THIRDS_PI*getRad()*getRad()*getRad()
            :
              M_PI*getRad()*getRad()
          );
      }
      
      void read(std::istream &iStream)
      {
        iStream >> m_pos >> m_rad;
      }

      void write(std::ostream &oStream) const
      {
        oStream
          << m_pos << " "
          << m_rad;
      }

    private:
      Vec3 m_pos;
      double m_rad;
      
      static bool s_is3d;
    };
  }
}

namespace std {
inline std::istream &operator>>(std::istream &iStream, esys::lsm::ParticleData &pd)
{
  pd.read(iStream);
  return iStream;
}

inline std::ostream &operator<<(std::ostream &oStream, const esys::lsm::ParticleData &pd)
{
  pd.write(oStream);
  return oStream;
}
}

namespace esys
{
  namespace lsm
  { 
    class Contact
    {
    public:
      Contact()
      {
      }

      Contact(
        const ParticleData &pd1,
        const ParticleData &pd2,
        const Vec3 &forcePos,
        const Vec3 &force
      )
        : m_pd1(pd1),
          m_pd2(pd2),
          m_forcePos(forcePos),
          m_force(force)
      {
      }

      Contact(const Contact &data)
        : m_pd1(data.getParticle1()),
          m_pd2(data.getParticle2()),
          m_forcePos(data.getForcePos()),
          m_force(data.getForce())
      {
      }

      const ParticleData &getParticle1() const
      {
        return m_pd1;
      }

      const ParticleData &getParticle2() const
      {
        return m_pd2;
      }

      const Vec3 &getCentrePos1() const
      {
        return getParticle1().getPos();
      }

      double getVolume1() const
      {
        return getParticle1().getVolume();
      }
            
      const Vec3 &getCentrePos2() const
      {
        return getParticle2().getPos();
      }

      double getVolume2() const
      {
        return getParticle2().getVolume();
      }

      const Vec3 &getForcePos() const
      {
        return m_forcePos;
      }

      const Vec3 &getForce() const
      {
        return m_force;
      }
          
      bool operator==(const Contact &data) const
      {
        return
          (
            (m_pd1 == data.m_pd1)
            &&
            (m_pd2 == data.m_pd2)
            &&
            (m_forcePos == data.m_forcePos)
            &&
            (m_force == data.m_force)
          );
      }

      void write(std::ostream &oStream) const
      {
        oStream
          << m_pd1      << " "
          << m_pd2      << " "
          << m_forcePos << " "
          << m_force;
      }

      void read(std::istream &iStream)
      {
        iStream
          >> m_pd1
          >> m_pd2
          >> m_forcePos
          >> m_force;
      }

    private:
      ParticleData m_pd1;
      ParticleData m_pd2;
      Vec3         m_forcePos;
      Vec3         m_force;
    };
  }
}

inline std::istream &operator>>(std::istream &iStream, esys::lsm::Contact &contact)
{
  contact.read(iStream);
  return iStream;
}

inline std::ostream &operator<<(std::ostream &oStream, const esys::lsm::Contact &contact)
{
  contact.write(oStream);
  return oStream;
}

#endif
