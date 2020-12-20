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


#ifndef ESYS_LSMCLOSEPACKITERATOR_H
#define ESYS_LSMCLOSEPACKITERATOR_H

#include "Foundation/BoundingBox.h"
#include "Foundation/vec3.h"
#include "Geometry/Vec3L.h"
#include "Geometry/ClosePackOrientation.h"

namespace esys
{
  namespace lsm
  {
    template <int NI, int NJ, int NK>
    class TmplMatrix
    {
    public:
      TmplMatrix();

      TmplMatrix(const TmplMatrix &m);

      TmplMatrix &operator=(const TmplMatrix &m);

      const double &operator()(int i, int j, int k) const;

      double &operator()(int i, int j, int k);

      int getNumI() const;

      int getNumJ() const;

      int getNumK() const;

    private:
      double m_matrix[NI][NJ][NK];
    };

    /**
     * Base class for iterators used to generate centre-points of spheres
     * arranged in a close-packing.
     */
    class ClosePackIterator
    {
    public:
      static const double SQRT_1_OVER_3;
      static const double SQRT_8_OVER_3;
      static const double SQRT_3;

      /**
       * Creates default empty iterator.
       */
      inline ClosePackIterator();

      /**
       * Creates an iterator which will iterate over numI*numJ*numK
       * centre points of spheres with radius sphereRadius.
       * @param numI number of spheres in the i direction.
       * @param numJ number of spheres in the j direction.
       * @param numK number of spheres in the k direction.
       * @param sphereRadius radius of spheres in the packing.
       * @param orientation specifies the axis alignment of layers.
       */
      inline ClosePackIterator(
        int numI,
        int numJ,
        int numK,
        double sphereRadius,
        ClosePackOrientation orientation = DEFAULT_ORIENT
      );

      /**
       * Returns whether there is another centre point in the
       * iteration sequence.
       */
      inline bool hasNext() const;

      /**
       * Returns the next centre-point in the iteration sequence.
       */
      inline Vec3 next();

      /**
       * Returns the radius of spheres used in the iteration.
       */
      inline double getRadius() const;

    protected:
      typedef TmplMatrix<3,6,6> OffsetMatrix;

      inline void incrementDimIndex();

      inline double getOffset(int i) const;

      inline const Vec3 &getMinPt() const;

      inline void setMinPt(const Vec3 &pt) const;

      inline void setDimRepeat(const Vec3L &dimRepeat);

      inline void setOffsetMatrix(const OffsetMatrix &offsetMatrix);

    private:
      static  Vec3L s_orientationDimMap[NUM_ORIENTATIONS];
      double        m_radius;
      Vec3          m_minPt;
      OffsetMatrix  m_offsetMatrix;
      Vec3L         m_dimRepeat;
      Vec3L         m_dimCount;
      Vec3L         m_dimIdx;
      Vec3L         m_dim;
    };
  }
}

#include "Geometry/ClosePackIterator.hpp"

#endif
