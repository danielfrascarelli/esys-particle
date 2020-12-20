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

#ifndef ESYS_LSMSIMPLEPARTICLEDATA_H
#define ESYS_LSMSIMPLEPARTICLEDATA_H

#include "Foundation/vec3.h"

#include <iostream>

namespace esys
{
  namespace lsm
  {
    /**
     * Container class for particle Id, Tag, position, radius and mass.
     */
    class SimpleParticleData
    {
    public:
      typedef int Id;
      typedef int Tag;

      inline SimpleParticleData();

      inline SimpleParticleData(Id id, Tag tag, const Vec3 &position, double radius);

      inline SimpleParticleData(const Vec3 &position, double radius, Id id, Tag tag);

      inline SimpleParticleData(const SimpleParticleData &p);

      inline SimpleParticleData &operator=(const SimpleParticleData &p);

      inline bool operator==(const SimpleParticleData &particleData) const;

      inline Id getId() const;

      inline void setId(const Id &id);

      inline Id getID() const;

      inline void setID(const Id &id);

      inline const Vec3 &getPosition() const;

      inline void setPosition(const Vec3 &pos);

      inline Tag getTag() const;

      inline void setTag(const Tag &tag);

      inline double getRadius() const;

      inline void setRadius(const double &r);

      inline void setMass(double mass);

      inline double getMass() const;

      inline double get2dMass() const;

      inline double get3dMass() const;

      inline void read(std::istream &istream);

      inline void write(std::ostream &write) const;

    private:
      Id     m_id;
      Tag    m_tag;
      Vec3   m_position;
      double m_radius;
      double m_mass;
    };
    inline std::istream &operator>>(std::istream &iStream, SimpleParticleData &particleData);
    inline std::ostream &operator<<(std::ostream &oStream, const SimpleParticleData &particleData);

  }
}

#include "Geometry/SimpleParticleData.hpp"

#endif
