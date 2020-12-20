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

namespace esys
{
  namespace lsm
  {
//===============================================================================
    SimpleParticleData::SimpleParticleData()
      : m_id(-1),
        m_tag(-1),
        m_position(),
        m_radius(0.0),
        m_mass(0.0)
    {
      setMass(get3dMass());
    }

    SimpleParticleData::SimpleParticleData(
      Id id,
      Tag tag,
      const Vec3 &position,
      double radius
    )
      : m_id(id),
        m_tag(tag),
        m_position(position),
        m_radius(radius),
        m_mass(0.0)
    {
      setMass(get3dMass());
    }

    SimpleParticleData::SimpleParticleData(
      const Vec3 &position,
      double radius,
      Id id,
      Tag tag
    )
      : m_id(id),
        m_tag(tag),
        m_position(position),
        m_radius(radius),
        m_mass(0.0)
    {
      setMass(get3dMass());
    }

    SimpleParticleData::SimpleParticleData(const SimpleParticleData &p)
      : m_id(p.m_id),
        m_tag(p.m_tag),
        m_position(p.m_position),
        m_radius(p.m_radius),
        m_mass(p.m_mass)
    {
    }

    SimpleParticleData &SimpleParticleData::operator=(const SimpleParticleData &p)
    {
      m_id =       p.m_id;
      m_tag =      p.m_tag;
      m_position = p.m_position;
      m_radius =   p.m_radius;
      m_mass =     p.m_mass;

      return *this;
    }

    bool SimpleParticleData::operator==(
      const SimpleParticleData &particleData
    ) const
    {
      return
        (
          (getId()       == particleData.getId())
          &&
          (getPosition() == particleData.getPosition())
          &&
          (getRadius()   == particleData.getRadius())
          &&
          (getTag()      == particleData.getTag())
        );
    }

    SimpleParticleData::Id SimpleParticleData::getId() const
    {
      return m_id;
    }

    void SimpleParticleData::setId(const Id &id)
    {
      m_id = id;
    }

    void SimpleParticleData::setID(const Id &id)
    {
      setId(id);
    }

    SimpleParticleData::Id SimpleParticleData::getID() const
    {
      return getId();
    }
    
    const Vec3 &SimpleParticleData::getPosition() const
    {
      return m_position;
    }

    void SimpleParticleData::setPosition(const Vec3 &pos)
    {
      m_position = pos;
    }

    SimpleParticleData::Tag SimpleParticleData::getTag() const
    {
      return m_tag;
    }

    void SimpleParticleData::setTag(const SimpleParticleData::Tag &tag)
    {
      m_tag = tag;
    }

    double SimpleParticleData::getRadius() const
    {
      return m_radius;
    }

    void SimpleParticleData::setRadius(const double &r)
    {
      m_radius = r;
    }

    void SimpleParticleData::setMass(double mass)
    {
      m_mass = mass;
    }
    
    double SimpleParticleData::getMass() const
    {
      return m_mass;
    }

    double SimpleParticleData::get2dMass() const
    {
      return m_radius*m_radius;
    }

    double SimpleParticleData::get3dMass() const
    {
      return m_radius*m_radius*m_radius;
    }

    void SimpleParticleData::read(std::istream &istream)
    {
      istream
        >> m_position
        >> m_radius
        >> m_id
        >> m_tag;
    }

    void SimpleParticleData::write(std::ostream &oStream) const
    {
      const char delim = ' ';
      oStream
        << getPosition() << delim
        << getRadius()   << delim
        << getId()       << delim
        << getTag();
    }

    std::istream &operator>>(std::istream &iStream, SimpleParticleData &particleData)
    {
      particleData.read(iStream);
      return iStream;
    }

    std::ostream &operator<<(std::ostream &oStream, const SimpleParticleData &particleData)
    {
      particleData.write(oStream);
      return oStream;
    }
  }
}
