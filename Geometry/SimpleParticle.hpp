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


SimpleParticle::SimpleParticle(const Vec3& pos,double rad,int id, int tag)
  : SimpleParticleData(pos, rad, id, tag)
{
}

SimpleParticle::SimpleParticle(const SimpleParticle& p)
  : SimpleParticleData(p)
{
}

SimpleParticle &SimpleParticle::operator=(const SimpleParticle& p)
{
  SimpleParticleData::operator=(p);
  return *this;
}

const Vec3 &SimpleParticle::getPos() const
{
  return getPosition();
}

void SimpleParticle::setPos(const Vec3 &pos)
{
  setPosition(pos);
}

void SimpleParticle::moveTo(const Vec3 &v)
{
  setPosition(v);
}

void SimpleParticle::translateBy(const Vec3 &v)
{
  setPosition(getPosition()+v);
}

void SimpleParticle::moveBy(const Vec3 &v)
{
  translateBy(v);
}

void SimpleParticle::rotate(const Vec3 &rotation, const Vec3 &posn)
{
  // From http://mathworld.wolfram.com/RotationFormula.html
  const double phi = rotation.norm();
  if (phi > 0.0)
  {
    const Vec3 r = getPosition() - posn;
    const Vec3 n = rotation/phi;
    const double cosPhi = cos(phi);
    const Vec3 rotatedR =
      r*cosPhi + n*((dot(n, r))*(1-cosPhi)) + cross(r, n)*sin(phi);
    setPosition(rotatedR + posn);
  }
}

double SimpleParticle::getRad() const
{
  return getRadius();
}

void SimpleParticle::setRad(double r)
{
  setRadius(r);
}

bool SimpleParticle::isValid() const
{
  return (getID() >= 0);
}

template <typename TmplVisitor>
void SimpleParticle::visit(const TmplVisitor &visitor) const
{
  visitor.visitSimpleParticle(*this);
}

template <typename TmplVisitor>
void SimpleParticle::visit(TmplVisitor &visitor)
{
  visitor.visitSimpleParticle(*this);
}

ostream& operator<<(ostream& ost,const SimpleParticle& p)
{
  ost
    << "Particle- id " << p.getId()
    << " pos: " << p.getPosition()
    << " rad: " << p.getRadius()
    << " tag : " << p.getTag() << std::endl;
  return ost;
}

ParticleComparer::ParticleComparer(const SimpleParticle &particle) : m_pParticle(&particle)
{
}

/*!
  \param p1
  \param p2
*/
bool ParticleComparer::operator()(const SimpleParticle &p1, const SimpleParticle &p2) const
{
  return (((p1.getPos() - m_pParticle->getPos()).norm() - p1.getRad())<
            ((p2.getPos() - m_pParticle->getPos()).norm() - p2.getRad()));
}

/*!
  \param p1
  \param p2
*/ 
bool ParticleComparer::operator()(const SimpleParticle *p1, const SimpleParticle *p2) const
{
  return (*this)(*p1, *p2);
}
