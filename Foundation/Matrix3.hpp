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

#ifndef __MATRIX3_HPP
#define __MATRIX3_HPP

#include "Foundation/Matrix3.h"
#include "Foundation/vec3.h"


MATRIX3_INLINE Matrix3::Matrix3()
{
  for(int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
      m[i][j]=0;
    }
  }
}

MATRIX3_INLINE Matrix3::Matrix3(const Vec3& V1, const Vec3& V2,const Vec3& V3)
{
  m[0][0]=V1.data[0];
  m[0][1]=V2.data[0];
  m[0][2]=V3.data[0];
  m[1][0]=V1.data[1];
  m[1][1]=V2.data[1];
  m[1][2]=V3.data[1];
  m[2][0]=V1.data[2];
  m[2][1]=V2.data[2];
  m[2][2]=V3.data[2];
  
}

MATRIX3_INLINE Matrix3::Matrix3(const double a[3][3])
{
  m[0][0]=a[0][0];
  m[0][1]=a[0][1];
  m[0][2]=a[0][2];
  m[1][0]=a[1][0];
  m[1][1]=a[1][1];
  m[1][2]=a[1][2];
  m[2][0]=a[2][0];
  m[2][1]=a[2][1];
  m[2][2]=a[2][2];
}

MATRIX3_INLINE Matrix3::Matrix3(const Matrix3& rhs)
{
  for(int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
      m[i][j]=rhs.m[i][j];
    }
  }
}

MATRIX3_INLINE Matrix3::~Matrix3()
{
}

/*!
  Determinant
*/
MATRIX3_INLINE double Matrix3::det()
{
  return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])+
    m[0][1]*(m[1][2]*m[2][0]-m[1][0]*m[2][2])+
    m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
}


MATRIX3_INLINE Matrix3 Matrix3::inv()
{
  Matrix3 res=*this;

  res.invert();
  
  return res;
}

MATRIX3_INLINE void Matrix3::transpose()
{
  double h;

  h=m[1][0];
  m[1][0]=m[0][1];
  m[0][1]=h;
  h=m[2][0];
  m[2][0]=m[0][2];
  m[0][2]=h;
  h=m[2][1];
  m[2][1]=m[1][2];
  m[1][2]=h;
}

MATRIX3_INLINE Matrix3 Matrix3::trans() const
{
  Matrix3 res;
  
  res.m[0][0]=m[0][0];
  res.m[0][1]=m[1][0];
  res.m[0][2]=m[2][0];
  res.m[1][0]=m[0][1];
  res.m[1][1]=m[1][1];
  res.m[1][2]=m[2][1];
  res.m[2][0]=m[0][2];
  res.m[2][1]=m[1][2];
  res.m[2][2]=m[2][2];

  return res;
}

// 15 Flops ( 9 mult,6 add)
MATRIX3_INLINE Vec3 Matrix3::operator *(const Vec3& V) const
{
  double x=m[0][0]*V.data[0]+m[0][1]*V.data[1]+m[0][2]*V.data[2];
  double y=m[1][0]*V.data[0]+m[1][1]*V.data[1]+m[1][2]*V.data[2];
  double z=m[2][0]*V.data[0]+m[2][1]*V.data[1]+m[2][2]*V.data[2];
  
  return Vec3(x,y,z);
}

// 9 Flops
MATRIX3_INLINE Matrix3 Matrix3::operator *(double d) const
{
  Matrix3 res;

  res.m[0][0]=m[0][0]*d;
  res.m[0][1]=m[0][1]*d;
  res.m[0][2]=m[0][2]*d;
  res.m[1][0]=m[1][0]*d;
  res.m[1][1]=m[1][1]*d;
  res.m[1][2]=m[1][2]*d;
  res.m[2][0]=m[2][0]*d;
  res.m[2][1]=m[2][1]*d;
  res.m[2][2]=m[2][2]*d;

  return res;
}

// 9 Flops
MATRIX3_INLINE Matrix3 Matrix3::operator /(double d) const
{
  Matrix3 res;

  res.m[0][0]=m[0][0]/d;
  res.m[0][1]=m[0][1]/d;
  res.m[0][2]=m[0][2]/d;
  res.m[1][0]=m[1][0]/d;
  res.m[1][1]=m[1][1]/d;
  res.m[1][2]=m[1][2]/d;
  res.m[2][0]=m[2][0]/d;
  res.m[2][1]=m[2][1]/d;
  res.m[2][2]=m[2][2]/d;

  return res;
}

MATRIX3_INLINE Matrix3 &Matrix3::operator=(const Matrix3 &other)
{
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      m[i][j] = other.m[i][j];
    }
  }
  return *this;
}

MATRIX3_INLINE bool Matrix3::operator==(const Matrix3 &other) const
{
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      if (m[i][j] != other.m[i][j]) {
        return false;
      }
    }
  }
  return true;
}

MATRIX3_INLINE std::ostream &operator<<(std::ostream &oStream, const Matrix3 &m)
{
  oStream << m(0,0);
  for(int i = 1; i < 9; i++) {
    oStream << " " << m(i/3, i%3);
  }
  return oStream;
}

// 45 Flops (27mult, 18add)
MATRIX3_INLINE Matrix3 Matrix3::operator *(const Matrix3& rhs) const
{
  Matrix3 res;

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      res.m[i][j]=m[i][0]*rhs.m[0][j]+m[i][1]*rhs.m[1][j]+m[i][2]*rhs.m[2][j];
    }
  }
  return res;
}

// 9 Flops
MATRIX3_INLINE Matrix3& Matrix3::operator+=(const Matrix3& rhs)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      m[i][j]=m[i][j]+rhs.m[i][j];
    }
  }
  return *this;
}

/*!
  add two matrices
*/
MATRIX3_INLINE Matrix3 Matrix3::operator +(const Matrix3& rhs) const
{
  Matrix3 res;
  
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      res.m[i][j]=m[i][j]+rhs.m[i][j];
    }
  }

  return res;
}

/*!
  subtract two matrices
*/
MATRIX3_INLINE Matrix3 Matrix3::operator-(const Matrix3& rhs) const
{
  Matrix3 res;
  
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) {
      res.m[i][j]=m[i][j]-rhs.m[i][j];
    }
  }

  return res;
}

/*!
  calculate the trace of a matrix
*/
MATRIX3_INLINE double Matrix3::trace() const
{
  return m[0][0]+m[1][1]+m[2][2];
}

/*!
  calculate the euclidian norm of a matrix
*/ 
MATRIX3_INLINE double Matrix3::norm() const
{
  double res=0.0;

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      res+=m[i][j]*m[i][j];
    }
  }
  return res;
}

// matrix from vector
MATRIX3_INLINE Matrix3 star(const Vec3& V)
{
  Matrix3 res;

  res.m[0][1]=-V.Z();
  res.m[0][2]=V.Y();
  res.m[1][0]=V.Z();
  res.m[1][2]=-V.X();
  res.m[2][0]=-V.Y();
  res.m[2][1]=V.X();

  return res;
}

// unit matrix
MATRIX3_INLINE Matrix3 Matrix3::Unit()
{
  Matrix3 res;

  res.m[0][0]=1.0;
  res.m[1][1]=1.0;
  res.m[2][2]=1.0;

  return res;
}

/*!
  scalar * matrix 
*/ 
MATRIX3_INLINE Matrix3 operator*(double d,const Matrix3& M)
{
  Matrix3 res;

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      res.m[i][j]=d*M.m[i][j];
    }
  }

  return res;
}



#endif // __MATRIX3_HPP
