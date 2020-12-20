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

#ifndef _QUATERNION_H
#define _QUATERNION_H

#define DO_INLINE_QUATERNION 1

#if DO_INLINE_QUATERNION >= 1
#define QUATERNION_INLINE inline
#else
#define QUATERNION_INLINE
#endif

#include <math.h>
#include "Foundation/vec3.h"

class Matrix3;

class Quaternion
{
private:
  Vec3   vector;
  double scalar;

public:
  //  Constructors
  QUATERNION_INLINE Quaternion();
  QUATERNION_INLINE Quaternion(double, const Vec3 &);

  //  Copy Constructor
  QUATERNION_INLINE Quaternion(const Quaternion &);

  //  Destructor
  QUATERNION_INLINE ~Quaternion() {};

  //  Assignment
  QUATERNION_INLINE Quaternion& operator=(const Quaternion&);

  //  Output
  QUATERNION_INLINE std::ostream& output(std::ostream&) const;
  QUATERNION_INLINE std::istream& input(std::istream& ci);

  //  Math
  QUATERNION_INLINE bool operator==(const Quaternion&) const;
  QUATERNION_INLINE bool operator!=(const Quaternion&) const;

  QUATERNION_INLINE Quaternion operator+(const Quaternion&) const;
  QUATERNION_INLINE Quaternion operator-(const Quaternion&) const;
  QUATERNION_INLINE Quaternion operator-() const;
  QUATERNION_INLINE friend Quaternion operator*(double, const Quaternion&);
  QUATERNION_INLINE Quaternion operator*(double) const;
  QUATERNION_INLINE Quaternion operator*(const Quaternion&) const;
  QUATERNION_INLINE Quaternion operator/(const Quaternion&) const;

  QUATERNION_INLINE Quaternion& operator+=(const Quaternion&);
  QUATERNION_INLINE Quaternion& operator-=(const Quaternion&);
  QUATERNION_INLINE Quaternion& operator*=(double);
  QUATERNION_INLINE Quaternion& operator*=(const Quaternion&);
  QUATERNION_INLINE Quaternion& operator/=(const Quaternion&);

  QUATERNION_INLINE Quaternion inverse() const;

  QUATERNION_INLINE void normalize();

  QUATERNION_INLINE double length() const;

  QUATERNION_INLINE Matrix3 to_matrix() const;

  //  Access Functions
  QUATERNION_INLINE Vec3 return_vec() const { return vector; };
  QUATERNION_INLINE double return_sca() const { return scalar; };

  QUATERNION_INLINE void set_vector(const Vec3 &v)  { vector = v; }
  QUATERNION_INLINE void set_scalar(double d) { scalar = d; }

  /**
   * Returns the angle and axis of rotation associated with this
   * quaternion as 3x1 vector. The magnitude of the vector is the
   * angle of rotation in radians.
   */
  QUATERNION_INLINE Vec3 asAngleAxis() const;

  /**
   * Pair representing angle of rotation about an axis.
   */
  typedef std::pair<double,Vec3> AngleAxisPair;
  /**
   * Returns the angle and axis of rotation associated with this
   * quaternion as std::pair<radians,3x1 vector>. Axis has non-unit
   * magnitude.
   */
  QUATERNION_INLINE AngleAxisPair asAngleAxisPair() const;
};

QUATERNION_INLINE std::ostream& operator<<(std::ostream&, const Quaternion &);
QUATERNION_INLINE std::istream& operator>>(std::istream&, Quaternion &);

#if DO_INLINE_QUATERNION >= 1
#include "Foundation/Quaternion.hpp"
#endif

#endif
