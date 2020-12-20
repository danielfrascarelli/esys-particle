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


#ifndef ESYS_LSMGRIDITERATOR_H
#define ESYS_LSMGRIDITERATOR_H

#include "Foundation/BoundingBox.h"
#include "Foundation/vec3.h"

namespace esys
{
  namespace lsm
  {
    /**
     * Class for iterating over the centre-points of spheres arranged
     * in a regular lattice within a specified bounding box.
     */
    class GridIterator
    {
    public:
      inline GridIterator()
        : m_sphereRadius(0.0),
          m_minI(0),
          m_minJ(0),
          m_minK(0),
          m_maxI(0),
          m_maxJ(0),
          m_maxK(0),
          m_i(0),
          m_j(0),
          m_k(0),
          m_centrePtBBox()
      {
      }

      inline GridIterator(int numPtsX, int numPtsY, int numPtsZ, double
sphereRadius)
        : m_sphereRadius(sphereRadius),
          m_minI(0),
          m_minJ(0),
          m_minK(0),
          m_maxI(numPtsX),
          m_maxJ(numPtsZ),
          m_maxK(numPtsY),
          m_i(0),
          m_j(0),
          m_k(0),
          m_centrePtBBox()
      {
        Vec3 minPt;
        Vec3 maxPt;
        if (numPtsZ > 1)
        {
          minPt = Vec3(0,0,0) + sphereRadius;
          maxPt = 
            Vec3(
              (numPtsX-1)*2.0*sphereRadius
              +
              (
                (numPtsZ > 1)
                ?
                sphereRadius
                :
                0.0
              )
              +
              (
                (numPtsY > 1)
                ?
                sphereRadius
                :
                0.0
              ),
              (numPtsY-1)*(sphereRadius*2.0*sqrt(2.0/3.0)),
              (numPtsZ-1)*(sphereRadius*sqrt(3.0))
              +
              ((numPtsY > 1) ? sphereRadius*sqrt(3.0)/3.0 : 0.0)
            );
        } else {
          minPt = Vec3(sphereRadius, sphereRadius, 0);
          maxPt = 
            Vec3(
              (numPtsX-1)*2.0*sphereRadius
              +
              (
                (numPtsY > 1)
                ?
                sphereRadius :
                0.0
              ),
              (numPtsY-1)*(sphereRadius*sqrt(3.0)),
              0
            );
          m_minJ = m_maxJ;
        }
        m_centrePtBBox = BoundingBox(minPt, (minPt + maxPt));
      }

	/*!
	  setup GridIterator from bounding box and sphere radius
	*/
	inline GridIterator(const BoundingBox &particleBBox, double sphereRadius, bool hard_limit=false)
        : m_sphereRadius(sphereRadius),
          m_minI(0),
          m_minJ(0),
          m_minK(0),
          m_maxI(0),
          m_maxJ(0),
          m_maxK(0),
          m_i(0),
          m_j(0),
          m_k(0),
          m_centrePtBBox()
      {
	int numPtsX ,numPtsY ,numPtsZ;
	if(!hard_limit){ // find "best fit"
	  numPtsX = max(1, int(nearbyint((particleBBox.getSizes().X()-(sphereRadius/4.0))/(sphereRadius*2.0))));
	  numPtsY = max(1, int(nearbyint(particleBBox.getSizes().Y()/(sphereRadius*2.0*sqrt(2.0/3.0)))));
	  numPtsZ = max(1, int(nearbyint(particleBBox.getSizes().Z()/(sphereRadius*sqrt(3.0)))));
	} else { // bounding box is hard limit
	  numPtsX = max(1, int(floor((particleBBox.getSizes().X()-(sphereRadius/4.0))/(sphereRadius*2.0))));
	  numPtsY = max(1, int(floor(particleBBox.getSizes().Y()/(sphereRadius*2.0*sqrt(2.0/3.0)))));
	  numPtsZ = max(1, int(floor(particleBBox.getSizes().Z()/(sphereRadius*sqrt(3.0)))));
	}
        if ((numPtsZ > 1) && (numPtsY > 1)) {
          numPtsX--;
        }
        if (particleBBox.getSizes().Z() > 0.0) {
          const Vec3 minPt = particleBBox.getMinPt() + sphereRadius;
          const Vec3 maxPt = 
            Vec3(
              (numPtsX-1)*2.0*sphereRadius + ((numPtsZ > 1) ? sphereRadius : 0.0) + ((numPtsY > 1) ? sphereRadius : 0.0),
              (numPtsY-1)*(sphereRadius*2.0*sqrt(2.0/3.0)),
              (numPtsZ-1)*(sphereRadius*sqrt(3.0))
              +
              ((numPtsY > 1) ? sphereRadius*sqrt(3.0)/3.0 : 0.0)
            );
          m_centrePtBBox = BoundingBox(minPt, (minPt + maxPt));
        } else {
          numPtsX = int(nearbyint((particleBBox.getSizes().X()-(sphereRadius/4.0))/(sphereRadius*2.0)));
          numPtsY = int(nearbyint(particleBBox.getSizes().Y()/(sphereRadius*sqrt(3.0))));
          numPtsZ = 0;

          Vec3 minPt = particleBBox.getMinPt() + sphereRadius;
          minPt.Z() = particleBBox.getMinPt().Z();
          const Vec3 maxPt = 
            Vec3(
              (numPtsX-1)*2.0*sphereRadius + ((numPtsY > 1) ? sphereRadius : 0.0),
              (numPtsY-1)*(sphereRadius*sqrt(3.0)),
              particleBBox.getMaxPt().Z()
            );
          m_centrePtBBox = BoundingBox(minPt, (minPt + maxPt));
        }
        m_minI = 0;
        m_minJ = 0;
        m_minK = 0;
        m_maxI = numPtsX;
        m_maxK = numPtsY;          
        m_maxJ = numPtsZ;

        m_i = m_minI;
        m_j = m_minJ;
        m_k = m_minK;
      }

      inline ~GridIterator()
      {
      }

      inline const BoundingBox &getBoundingBox() const
      {
        return m_centrePtBBox;
      }

      inline BoundingBox getSphereBBox() const
      {
        return
          (
            is2d()
            ?
              BoundingBox(
                m_centrePtBBox.getMinPt() - Vec3(m_sphereRadius, m_sphereRadius, 0.0),
                m_centrePtBBox.getMaxPt() + Vec3(m_sphereRadius, m_sphereRadius, 0.0)
              )
            :
              BoundingBox(
                m_centrePtBBox.getMinPt() - m_sphereRadius,
                m_centrePtBBox.getMaxPt() + m_sphereRadius
              )
          );
      }

      /**
       * Returns whether there are any points remaining
       * in the iteration.
       */
      inline bool hasNext() const
      {
        return (m_i < m_maxI);
      }

      inline bool is2d() const
      {
        return (m_minJ == m_maxJ);
      }

      inline Vec3 getPoint() const
      {
        if (is2d()) {
          return 
            Vec3(
              m_centrePtBBox.getMinPt().X() + (double(m_i)+0.5*double(m_k%2))*m_sphereRadius*2.0,
              m_centrePtBBox.getMinPt().Y() + double(m_k)*sqrt(3.0)*m_sphereRadius,
              0.0
            );
        } else {
          return
            Vec3(
              m_centrePtBBox.getMinPt().X() + (double(m_i)+0.5*double(m_j%2)+0.5*double(m_k%2))*m_sphereRadius*2.0,
              m_centrePtBBox.getMinPt().Y() + (double(m_k)*2.0*sqrt(2.0/3.0))*m_sphereRadius,
              m_centrePtBBox.getMinPt().Z() + ((double(m_j)+double(m_k%2)/3.0)*sqrt(3.0))*m_sphereRadius
            );
        }
      }

      inline void increment()
      {
        if (m_k+1 < m_maxK) {
          m_k++;
        }
        else if (m_j+1 < m_maxJ) {
          m_k = m_minK;
          m_j++;
        }
        else {
          m_k = m_minK;
          m_j = m_minJ;
          m_i++;
        }
      }

      /**
       * Returns the next regular-grid-point in the iteration.
       */
      inline Vec3 next()
      {
        const Vec3 pt = getPoint();
        increment();
        return pt;
      }

    private:
      double m_sphereRadius;

      int m_minI;
      int m_minJ;
      int m_minK;

      int m_maxI;
      int m_maxJ;
      int m_maxK;

      int m_i; //! index in x-direction
      int m_j; //! index in y-direction
      int m_k; //! index in z-direction
      
      BoundingBox m_centrePtBBox;
    };
  }
}

#endif
