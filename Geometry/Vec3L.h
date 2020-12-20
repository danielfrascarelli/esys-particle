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


#ifndef ESYS_LSMVEC3L_H
#define ESYS_LSMVEC3L_H

#include <iostream>

namespace esys
{

  namespace lsm
  {

    /**
    @author Shane J Latham
    */
    class Vec3L
    {
    public:
      typedef long Long;

      Vec3L()
      {
        m_data[0] = 0;
        m_data[1] = 0;
        m_data[2] = 0;
      }
      
      Vec3L(Long x, Long y, Long z)
      {
        m_data[0] = x;
        m_data[1] = y;
        m_data[2] = z;
      }
      
      Vec3L(const Vec3L &vec)
      {
        m_data[0] = vec.m_data[0];
        m_data[1] = vec.m_data[1];
        m_data[2] = vec.m_data[2];
      }

      ~Vec3L()
      {
      }

      Vec3L &operator=(const Vec3L &vec)
      {
        m_data[0] = vec.m_data[0];
        m_data[1] = vec.m_data[1];
        m_data[2] = vec.m_data[2];
        
        return *this;
      }

      bool operator==(const Vec3L &vec) const
      {
        return
          (
            (m_data[0] == vec.m_data[0])
            &&
            (m_data[1] == vec.m_data[1])
            &&
            (m_data[2] == vec.m_data[2])
          );
      }

      Long &operator[](int idx)
      {
        return m_data[idx];
      }

      const Long &operator[](int idx) const
      {
        return m_data[idx];
      }

      Vec3L operator-(Long val) const
      {
        return Vec3L(m_data[0]-val, m_data[1]-val, m_data[2]-val);
      }

      Vec3L operator+(Long val) const
      {
        return Vec3L(m_data[0]+val, m_data[1]+val, m_data[2]+val);
      }

      const Long &X() const
      {
        return m_data[0];
      }

      Long &X()
      {
        return m_data[0];
      }

      const Long &Y() const
      {
        return m_data[1];
      }

      Long &Y()
      {
        return m_data[1];
      }
      
      const Long &Z() const
      {
        return m_data[2];
      }
      
      Long &Z()
      {
        return m_data[2];
      }
      
      Vec3L min(const Vec3L &vec) const
      {
        return
          Vec3L
          (
            std::min(m_data[0], vec.m_data[0]),
            std::min(m_data[1], vec.m_data[1]),
            std::min(m_data[2], vec.m_data[2])
          );
      }

      Vec3L max(const Vec3L &vec) const
      {
        return
          Vec3L
          (
            std::max(m_data[0], vec.m_data[0]),
            std::max(m_data[1], vec.m_data[1]),
            std::max(m_data[2], vec.m_data[2])
          );
      }

    private:
      Long m_data[3];
    };

    inline std::ostream &operator<<(std::ostream &oStream, const Vec3L &vec)
    {
      oStream << vec.X() << " " << vec.Y() << " " << vec.Z();
      return oStream;
    }
  }
}

#endif
