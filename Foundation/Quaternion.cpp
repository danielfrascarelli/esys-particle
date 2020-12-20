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

//
// ============================================================
//
//                      Quaternion.cpp

#include "Foundation/Quaternion.h"
#include "Foundation/Quaternion.hpp"

#if 0
//
// ============================================================
//
//             CONSTRUCTORS, DESTRUCTORS
//
// ============================================================
//

Quaternion::Quaternion()
  : vector(Vec3::ZERO),
    scalar(1.0)
{
}

Quaternion::Quaternion(double d, const Vec3 &v)
  : vector(v),
    scalar(d)
{
}

Quaternion::Quaternion(const Quaternion & q)
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

Quaternion& Quaternion::operator=(const Quaternion& q)
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

std::ostream& operator<<(std::ostream& co, const Quaternion& q)
{
  return q.output(co);
}

std::istream& operator>>(std::istream& ci, Quaternion& q)
{
  return q.input(ci);
}

std::ostream& Quaternion::output(std::ostream& co) const 
{
  co
     << scalar << ' '
     << vector;

  return co;
}

std::istream& Quaternion::input(std::istream& ci)
{
  ci
     >> scalar
     >> vector;

  return ci;
}

bool Quaternion::operator==(const Quaternion &q) const
{
  return
    (
      (return_sca() == q.return_sca())
      &&
      (return_vec() == q.return_vec())
    );
}

bool Quaternion::operator!=(const Quaternion &q) const
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

Quaternion Quaternion::operator+(const Quaternion& q2) const
{
#if 0
  Quaternion qq;

  qq.scalar = scalar + q2.scalar;
  qq.vector = vector + q2.vector;
#endif
  return Quaternion(scalar + q2.scalar, vector + q2.vector);
}

Quaternion Quaternion::operator-(const Quaternion& q2) const
{
#if 0
  Quaternion qq;

  qq.scalar = scalar - q2.scalar;
  qq.vector = vector - q2.vector;
#endif
  return Quaternion(scalar-q2.scalar, vector-q2.vector);
}

Quaternion Quaternion::operator-() const
{
#if 0
  Quaternion qq;

  qq.scalar = -scalar;
  qq.vector = -vector;
#endif
  return Quaternion(-scalar, -vector);
}

Quaternion operator*(double c, const Quaternion& q) 
{
#if 0
  Quaternion qq;

  qq.scalar = c * q.scalar;
  qq.vector = c * q.vector;
#endif
  return Quaternion(c*q.scalar, c*q.vector);
}

Quaternion Quaternion::operator*(double c) const
{
  return Quaternion(c * scalar, c * vector);
}

Quaternion Quaternion::operator*(const Quaternion& q2) const
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

Quaternion Quaternion::inverse() const
{
  return Quaternion(scalar, -vector);
}

Quaternion Quaternion::operator/(const Quaternion& q2) const
{
  return *this * q2.inverse();
}

Quaternion& Quaternion::operator+=(const Quaternion& q)
{
  scalar += q.scalar;
  vector += q.vector;

  return *this;
}

Quaternion& Quaternion::operator-=(const Quaternion& q)
{
  scalar -= q.scalar;
  vector -= q.vector;

  return *this;
}

Quaternion& Quaternion::operator*=(double c)
{
  scalar *= c;
  vector *= c;

  return *this;
}

Quaternion& Quaternion::operator*=(const Quaternion& q)
{
  const double s = scalar * q.scalar - dot(vector, q.vector);
  vector = scalar * q.vector
          + q.scalar * vector
          + cross(vector, q.vector);
  scalar = s;

  return *this;
}

Quaternion& Quaternion::operator/=(const Quaternion& q)
{
  Quaternion qq = q.inverse();

  const double s = scalar * qq.scalar - dot(vector, qq.vector);
  vector = scalar * qq.vector
          + qq.scalar * vector
          + cross(vector, qq.vector);
  scalar = s;

  return *this;
}


void Quaternion::normalize() 
{
  double len = length();

  if (len > 1.0e-8)
  {
    scalar /= len;
    vector /= len;
  }
}

double Quaternion::length() const 
{
  const double vlen = vector.norm();
  return sqrt(scalar * scalar + vlen * vlen);
}


Matrix3 Quaternion::to_matrix() const
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

Vec3 Quaternion::asAngleAxis() const
{
  Vec3 v(this->vector);
  return (v *= (2*acos(this->scalar)/v.norm()));
}

Quaternion::AngleAxisPair Quaternion::asAngleAxisPair() const
{
  return AngleAxisPair(2*acos(this->scalar), this->vector);
}

#endif // 0
