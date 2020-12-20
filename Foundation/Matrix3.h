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

#ifndef __MATRIX3_H
#define __MATRIX3_H

#define DO_INLINE_MATRIX3 1

#if DO_INLINE_MATRIX3 >= 1
#define MATRIX3_INLINE inline
#else
#define MATRIX3_INLINE
#endif

// --- IO includes ---
#include <iostream>
using std::ostream;

// --- STL includes ---
#include <array>

//! exception class for Matrix3
class MatSingularError
{
public:
  MATRIX3_INLINE MatSingularError()
  {
  };
};


class Vec3 ;
/*!
  \class Matrix3
  \brief 3x3 Matrix
 
  \author Steffen Abe
  $Revision$
  $Date$ 
*/
class Matrix3
{
private:
  double m[3][3];

public:
  MATRIX3_INLINE Matrix3();
  MATRIX3_INLINE Matrix3(const Vec3&, const Vec3&,const Vec3&);
  MATRIX3_INLINE Matrix3(const double[3][3]);
  Matrix3(const double[9]);
  MATRIX3_INLINE Matrix3(const Matrix3&);
  MATRIX3_INLINE virtual ~Matrix3();

  MATRIX3_INLINE double det();
  Vec3 solve(const Vec3&) const;    
  Vec3 solve_homogeneous() const;
  void invert(); //!< in-situ inversion
  MATRIX3_INLINE Matrix3 inv(); //!< return inverse;
  MATRIX3_INLINE void transpose(); //!< transpose in situ
  MATRIX3_INLINE Matrix3 trans() const; //!< return transposed
  MATRIX3_INLINE Vec3 operator *(const Vec3&) const;
  MATRIX3_INLINE Matrix3 operator *(double) const;
  MATRIX3_INLINE Matrix3 operator /(double) const;
  MATRIX3_INLINE Matrix3 operator *(const Matrix3&) const;
  MATRIX3_INLINE Matrix3& operator +=(const Matrix3&);
  MATRIX3_INLINE Matrix3 operator +(const Matrix3&)const;
  MATRIX3_INLINE Matrix3 operator -(const Matrix3&)const;
  MATRIX3_INLINE bool operator==(const Matrix3&)const;
  MATRIX3_INLINE Matrix3 &operator=(const Matrix3&);
  MATRIX3_INLINE double trace() const;
  MATRIX3_INLINE double norm() const;
  Matrix3 get_symmetric_part() const;
  Matrix3 dot(const Matrix3&);

  MATRIX3_INLINE double operator()(int i, int j) const {return m[i][j];}

  MATRIX3_INLINE double& operator()(int i, int j){return m[i][j];}

  MATRIX3_INLINE friend Matrix3 operator*(double,const Matrix3&);

  //!< generate matrix from vector
  MATRIX3_INLINE friend Matrix3 star(const Vec3&);

  //!< generate unit matrix
  MATRIX3_INLINE static Matrix3 Unit();

  //!< eigenvectors, eigenvalues
  void eigen(Vec3&,Vec3&,Vec3&,double&,double&,double&);
  std::array<double,3> eigenvalues() const;

  // output
  MATRIX3_INLINE friend ostream& operator<<(ostream&,const Matrix3&);
}; 

#if DO_INLINE_MATRIX3 >= 1
#include "Foundation/Matrix3.hpp"
#endif

#endif //__MATRIX3_H
