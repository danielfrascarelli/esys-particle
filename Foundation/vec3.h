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

#ifndef __VEC3_H
#define __VEC3_H

#define DO_INLINE_VEC3 1

#if DO_INLINE_VEC3 >= 1
#define VEC3_INLINE inline
#else
#define VEC3_INLINE
#endif

#include <iostream>
#include <math.h>
#include <string>

#include "Foundation/Error.h"

using std::ostream;
using std::istream;
using std::string;

class Matrix3;

class VecErr:public MError
{
 public:
  VecErr(const string&);
  virtual ~VecErr(){};
};

struct VDMulVadd;
struct VDMul;

class Vec3
{
protected:
  double data[3];

public:
  static const Vec3 ZERO; //! The zero vector.
  // constructors
  VEC3_INLINE Vec3();
  VEC3_INLINE explicit Vec3(double s);
  VEC3_INLINE Vec3(double,double,double);
  VEC3_INLINE Vec3(const Vec3&);

  // vec-vec operators
  VEC3_INLINE Vec3& operator=(const Vec3&);
  VEC3_INLINE Vec3& operator=(double s);
  VEC3_INLINE Vec3& operator-=(const Vec3&);
  VEC3_INLINE Vec3& operator+=(const Vec3&);
  VEC3_INLINE Vec3 operator+(const Vec3&) const;
  VEC3_INLINE Vec3 operator-(const Vec3&) const;
//wangyc added !
  VEC3_INLINE Vec3 operator*(const Matrix3 &m) const;
  VEC3_INLINE double operator*(const Vec3&) const; 
  VEC3_INLINE Vec3 operator-() const;
  
  // vec-dbl ops
  VEC3_INLINE Vec3 operator*(double) const;
  VEC3_INLINE Vec3& operator*=(double);
  VEC3_INLINE Vec3 operator/(double) const;
  VEC3_INLINE Vec3 operator-(double) const;
  VEC3_INLINE Vec3 operator+(double) const;
  VEC3_INLINE Vec3& operator+=(double);
  VEC3_INLINE Vec3& operator-=(double);

// wangyc added !
  VEC3_INLINE Vec3& operator/=(double);
  VEC3_INLINE double norm() const;
  VEC3_INLINE double norm2() const;
  VEC3_INLINE Vec3 unit() const;
  VEC3_INLINE Vec3 unit_s() const; //safe version (throw exceptions)
  VEC3_INLINE double max() const;
  VEC3_INLINE double min() const;

  VEC3_INLINE Vec3 rotate(const Vec3 &axis, const Vec3 &axisPt) const;
  VEC3_INLINE void rotateBy(const Vec3& origin, const Vec3& axis, double angle); 

  VEC3_INLINE bool operator==(const Vec3&) const;
  VEC3_INLINE bool operator!=(const Vec3&) const;

  VEC3_INLINE friend Vec3 cmax(const Vec3&,const Vec3&);
  VEC3_INLINE friend Vec3 cmin(const Vec3&,const Vec3&);

  VEC3_INLINE friend Vec3 cross(const Vec3&,const Vec3&);
  VEC3_INLINE friend double dot(const Vec3&,const Vec3&); 
  VEC3_INLINE friend Vec3 operator*(double,const Vec3&);

  //n+1-ary operators
  VEC3_INLINE void mul_add_and_assign(const Vec3*,const Vec3*,const double&);
  VEC3_INLINE void mul_and_assign(const Vec3*,const double&);

  VEC3_INLINE Vec3(const VDMulVadd&);
  VEC3_INLINE Vec3& operator=(const VDMulVadd&);

  VEC3_INLINE Vec3(const VDMul&);
  VEC3_INLINE Vec3& operator=(const VDMul&);

  //access stuff
// wangyc added ! 
  VEC3_INLINE void set_x(double x) {data[0] = x;}
  VEC3_INLINE void set_y(double y) {data[1] = y;}
  VEC3_INLINE void set_z(double z) {data[2] = z;}
//  void set_xyz(double x, double y, double z)
//  { data[0] = x; data[1] = y; data[2] = z;}
  
  VEC3_INLINE double& X() {return data[0];};
  VEC3_INLINE double& Y() {return data[1];};
  VEC3_INLINE double& Z() {return data[2];};
  VEC3_INLINE double X() const {return data[0];};
  VEC3_INLINE double Y() const {return data[1];};
  VEC3_INLINE double Z() const {return data[2];};
  VEC3_INLINE const double &operator[](int i) const {return data[i];}
  VEC3_INLINE double& operator[](int i) {return data[i];}

  // in/output
  VEC3_INLINE friend ostream& operator << (ostream&,const Vec3&);
  VEC3_INLINE friend istream& operator >> (istream&,Vec3&);

  // comparison -> enable to use of Vec3 as key in STL map and set
  bool operator<(const Vec3&) const; 

  friend class Matrix3;
};

VEC3_INLINE Vec3 comp_max(const Vec3&,const Vec3&); //!< per component maximum
VEC3_INLINE Vec3 comp_min(const Vec3&,const Vec3&); //!< per component minimum

#if DO_INLINE_VEC3 >= 1
#include "Foundation/vec3.hpp"
#endif

#endif // __VEC3_H
