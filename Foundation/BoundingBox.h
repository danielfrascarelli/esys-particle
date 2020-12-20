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


#ifndef ESYS_LSMBOUNDINGBOX_H
#define ESYS_LSMBOUNDINGBOX_H

#include "Foundation/vec3.h"

namespace esys
{
  namespace lsm
  {
    /*!
     
      \brief 3D bounding box 
    */
    class BoundingBox
    {
    public:
      inline BoundingBox();
      
      inline BoundingBox(const Vec3 &minBBoxPt, const Vec3 &maxBBoxPt);

      inline virtual ~BoundingBox();

      inline double getVolume() const;

      inline const Vec3 &getMinPt() const;

      inline const Vec3 &getMaxPt() const;
      
      inline bool operator==(const BoundingBox &bbox) const;
      
      inline Vec3 getSizes() const;
      
      inline bool contains(const Vec3 &pt, double tolerance = 0.0) const;

    private:
      Vec3 m_minPt;
      Vec3 m_maxPt;
    };
    inline std::ostream &operator<<(std::ostream &oStream, const BoundingBox &bbox);
  }
}

#include "Foundation/BoundingBox.hpp"

#endif
