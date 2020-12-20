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


#ifndef _QUATERNION_HPP
#define _QUATERNION_HPP
//
// ============================================================
//
//                      Quaternion.hpp

#include "Foundation/Quaternion.h"
#include "Foundation/vec3.h"
#include "Foundation/Matrix3.h"

//
// ============================================================
//
//             CONSTRUCTORS, DESTRUCTORS
//
// ============================================================
//

QUATERNION_INLINE Quaternion::Quaternion()
  : vector(Vec3::ZERO),
    scalar(1.0)
{
}

QUATERNION_INLINE Quaternion::Quaternion(double d, const Vec3 &v)
  : vector(v),
    scalar(d)
{
}

QUATERNION_INLINE Quaternion::Quaternion(const Quaternion & q)
  : vector(q.vector),
    scalar(q.scalar)
{
}

//
// ============================================================
//
//            ASSIGNMENT
//
// ============================================================
//

QUATERNION_INLINE Quaternion& Quaternion::operator=(const Quaternion& q)
{
#if 0
  if (&q == this) return *this;
#endif
  vector = q.vector;
  scalar = q.scalar;

  return *this;
}

//
// ============================================================
//
//            OUTPUT
//
// ============================================================
//

QUATERNION_INLINE std::ostream& operator<<(std::ostream& co, const Quaternion& q)
{
  return q.output(co);
}

QUATERNION_INLINE std::istream& operator>>(std::istream& ci, Quaternion& q)
{
  return q.input(ci);
}

QUATERNION_INLINE std::ostream& Quaternion::output(std::ostream& co) const 
{
  co
     << scalar << ' '
     << vector;

  return co;
}

QUATERNION_INLINE std::istream& Quaternion::input(std::istream& ci)
{
  ci
     >> scalar
     >> vector;

  return ci;
}

QUATERNION_INLINE bool Quaternion::operator==(const Quaternion &q) const
{
  return
    (
      (return_sca() == q.return_sca())
      &&
      (return_vec() == q.return_vec())
    );
}

QUATERNION_INLINE bool Quaternion::operator!=(const Quaternion &q) const
{
  return !(*this == q);
}

//
// ============================================================
//
//            ARITHMETIC OPERATIONS
//
// ============================================================
//

QUATERNION_INLINE Quaternion Quaternion::operator+(const Quaternion& q2) const
{
#if 0
  Quaternion qq;

  qq.scalar = scalar + q2.scalar;
  qq.vector = vector + q2.vector;
#endif
  return Quaternion(scalar + q2.scalar, vector + q2.vector);
}

QUATERNION_INLINE Quaternion Quaternion::operator-(const Quaternion& q2) const
{
#if 0
  Quaternion qq;

  qq.scalar = scalar - q2.scalar;
  qq.vector = vector - q2.vector;
#endif
  return Quaternion(scalar-q2.scalar, vector-q2.vector);
}

QUATERNION_INLINE Quaternion Quaternion::operator-() const
{
#if 0
  Quaternion qq;

  qq.scalar = -scalar;
  qq.vector = -vector;
#endif
  return Quaternion(-scalar, -vector);
}

QUATERNION_INLINE Quaternion operator*(double c, const Quaternion& q) 
{
#if 0
  Quaternion qq;

  qq.scalar = c * q.scalar;
  qq.vector = c * q.vector;
#endif
  return Quaternion(c*q.scalar, c*q.vector);
}

QUATERNION_INLINE Quaternion Quaternion::operator*(double c) const
{
  return Quaternion(c * scalar, c * vector);
}

QUATERNION_INLINE Quaternion Quaternion::operator*(const Quaternion& q2) const
{
#if 0
  Quaternion qq;

  qq.scalar = scalar * q2.scalar - dot(vector, q2.vector);
  qq.vector = scalar * q2.vector
             + q2.scalar * vector
             + cross(vector, q2.vector);
#endif
  return 
    Quaternion(
      scalar * q2.scalar - dot(vector, q2.vector),
      scalar * q2.vector
      + q2.scalar * vector
      + cross(vector, q2.vector)
    );
}

QUATERNION_INLINE Quaternion Quaternion::inverse() const
{
  return Quaternion(scalar, -vector);
}

QUATERNION_INLINE Quaternion Quaternion::operator/(const Quaternion& q2) const
{
  return *this * q2.inverse();
}

QUATERNION_INLINE Quaternion& Quaternion::operator+=(const Quaternion& q)
{
  scalar += q.scalar;
  vector += q.vector;

  return *this;
}

QUATERNION_INLINE Quaternion& Quaternion::operator-=(const Quaternion& q)
{
  scalar -= q.scalar;
  vector -= q.vector;

  return *this;
}

QUATERNION_INLINE Quaternion& Quaternion::operator*=(double c)
{
  scalar *= c;
  vector *= c;

  return *this;
}

QUATERNION_INLINE Quaternion& Quaternion::operator*=(const Quaternion& q)
{
  const double s = scalar * q.scalar - dot(vector, q.vector);
  vector = scalar * q.vector
          + q.scalar * vector
          + cross(vector, q.vector);
  scalar = s;

  return *this;
}

QUATERNION_INLINE Quaternion& Quaternion::operator/=(const Quaternion& q)
{
  Quaternion qq = q.inverse();

  const double s = scalar * qq.scalar - dot(vector, qq.vector);
  vector = scalar * qq.vector
          + qq.scalar * vector
          + cross(vector, qq.vector);
  scalar = s;

  return *this;
}


QUATERNION_INLINE void Quaternion::normalize() 
{
  double len = length();

  if (len > 1.0e-8)
  {
    scalar /= len;
    vector /= len;
  }
}

QUATERNION_INLINE double Quaternion::length() const 
{
  const double vlen = vector.norm();
  return sqrt(scalar * scalar + vlen * vlen);
}

QUATERNION_INLINE Matrix3 Quaternion::to_matrix() const
{
 double m[3][3];

  m[0][0] =   scalar*scalar + vector.X()*vector.X()
            - vector.Y()*vector.Y() -vector.Z()*vector.Z();
  m[0][1] =   2*(vector.X()*vector.Y() + scalar*vector.Z());
  m[0][2] =   2*(vector.X()*vector.Z() - scalar*vector.Y());
  m[1][0] =   2*(vector.X()*vector.Y() - scalar*vector.Z());
  m[1][1] =   scalar*scalar - vector.X()*vector.X()
            + vector.Y()*vector.Y() - vector.Z()*vector.Z();
  m[1][2] =   2*(vector.Y()*vector.Z() + scalar*vector.X());
  m[2][0] =   2*(vector.X()*vector.Z() + scalar*vector.Y());
  m[2][1] =   2*(vector.Y()*vector.Z() - scalar*vector.X());
  m[2][2] =   scalar*scalar - vector.X()*vector.X()
            - vector.Y()*vector.Y() + vector.Z()*vector.Z();
 // In 2D case m[2][2] is not so accurate!
//    m[2][2] = 0.0 ;

/*

  m[0][0] =  1.0 
            - 2.0*vector.Y()*vector.Y() -2.0*vector.Z()*vector.Z();
  m[0][1] =   2*(vector.X()*vector.Y() + scalar*vector.Z());
  m[0][2] =   2*(vector.X()*vector.Z() - scalar*vector.Y());
  m[1][0] =   2*(vector.X()*vector.Y() - scalar*vector.Z());
  m[1][1] =  1.0  -2.0* vector.X()*vector.X()
             - 2.0*vector.Z()*vector.Z();
  m[1][2] =   2*(vector.Y()*vector.Z() + scalar*vector.X());
  m[2][0] =   2*(vector.X()*vector.Z() + scalar*vector.Y());
  m[2][1] =   2*(vector.Y()*vector.Z() - scalar*vector.X());
  m[2][2] =   1.0 - 2.0*vector.X()*vector.X()
            - 2.0*vector.Y()*vector.Y() ;
*/

  return  Matrix3(m);
}

QUATERNION_INLINE Vec3 Quaternion::asAngleAxis() const
{
  Vec3 v(this->vector);
  return (v *= (2*acos(this->scalar)/v.norm()));
}

QUATERNION_INLINE Quaternion::AngleAxisPair Quaternion::asAngleAxisPair() const
{
  return AngleAxisPair(2*acos(this->scalar), this->vector);
}

#endif // _QUATERNION_HPP
