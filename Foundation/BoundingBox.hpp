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
    BoundingBox::BoundingBox()
      : m_minPt(Vec3::ZERO),
        m_maxPt(Vec3::ZERO)
    {
    }

    BoundingBox::BoundingBox(const Vec3 &minBBoxPt, const Vec3 &maxBBoxPt)
      : m_minPt(minBBoxPt),
        m_maxPt(maxBBoxPt)
    {
    }

    BoundingBox::~BoundingBox()
    {
    }

    double BoundingBox::getVolume() const
    {
      Vec3 diff(m_maxPt-m_minPt);
      return fabs(diff.X()*diff.Y()*diff.Z());
    }

    const Vec3 &BoundingBox::getMinPt() const
    {
      return m_minPt;
    }

    const Vec3 &BoundingBox::getMaxPt() const
    {
      return m_maxPt;
    }

    bool BoundingBox::operator==(const BoundingBox &bbox) const
    {
      return
        (
          (getMinPt() == bbox.getMinPt())
          &&
          (getMaxPt() == bbox.getMaxPt())
        );
    }

    bool BoundingBox::contains(const Vec3 &pt, double tolerance) const
    {
      return
        (
          (getMinPt().X() <= (pt.X() + tolerance))
          &&
          (getMinPt().Y() <= (pt.Y() + tolerance))
          &&
          (getMinPt().Z() <= (pt.Z() + tolerance))
          &&
          (getMaxPt().X() >= (pt.X() - tolerance))
          &&
          (getMaxPt().Y() >= (pt.Y() - tolerance))
          &&
          (getMaxPt().Z() >= (pt.Z() - tolerance))
        );
    }
    
    Vec3 BoundingBox::getSizes() const
    {
      return (getMaxPt() - getMinPt());
    }

    std::ostream &operator<<(std::ostream &oStream, const BoundingBox &bbox)
    {
      oStream << bbox.getMinPt() << " " << bbox.getMaxPt();
      return oStream;
    }
  }
}
