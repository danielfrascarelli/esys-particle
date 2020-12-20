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

#include <iostream>
#include <math.h>
#include <string>

#include "Error.h"

using std::ostream;
using std::istream;
using std::string;

class Mat3;

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
  // constructors
  Vec3();
  Vec3(double,double,double);
  Vec3(const Vec3&);

  // vec-vec operators
  Vec3& operator=(const Vec3&);
  Vec3& operator-=(const Vec3&);
  Vec3& operator+=(const Vec3&);
  Vec3 operator+(const Vec3&) const;
  Vec3 operator-(const Vec3&) const;
  double operator*(const Vec3&) const; 
  inline Vec3 operator-() { return Vec3(-data[0],-data[1],-data[2]); } ; 
  
  // vec-dbl ops
  Vec3 operator*(double) const;
  Vec3 operator/(double) const;

  double norm() const;
  double norm2() const;
  Vec3 unit() const;
  Vec3 unit_s() const; //safe version (throw exceptions)
  double max() const;
  double min() const;

  bool operator==(const Vec3&);
  bool operator!=(const Vec3&);

  friend Vec3 cmax(const Vec3&,const Vec3&);
  friend Vec3 cmin(const Vec3&,const Vec3&);

  friend Vec3 cross(const Vec3&,const Vec3&);
  
  friend Vec3 operator*(double,const Vec3&);

  //n+1-ary operators
  void mul_add_and_assign(const Vec3*,const Vec3*,const double&);
  void mul_and_assign(const Vec3*,const double&);

  Vec3(const VDMulVadd&);
  Vec3& operator=(const VDMulVadd&);

  Vec3(const VDMul&);
  Vec3& operator=(const VDMul&);

  //access stuff
  inline double X() const {return data[0];};
  inline double Y() const {return data[1];};
  inline double Z() const {return data[2];};
  inline double operator[](int i) const {return data[i];}; 
  inline double& operator[](int i) {return data[i];}; 

  // in/output
  friend ostream& operator << (ostream&,const Vec3&);
  friend istream& operator >> (istream&,Vec3&);

  friend class Mat3;
};

//------------------------------
// stuff for 'n+1-ary' ops
// see Stroustrup, p.675 f.
//------------------------------

struct VDMul
{
  const Vec3& v;
  const double d;

  VDMul(const Vec3& vv, const double& dd):v(vv),d(dd){}

  //  operator Vec3(){return Vec3(v.x*d,v.y*d,v.z*d);};
};

inline VDMul operator*(const Vec3& vv, const double& dd)
{
  return VDMul(vv,dd);
}

struct VDMulVadd
{
  const Vec3& v1;
  const Vec3& v2;
  const double& d;

  VDMulVadd(const VDMul vd,const Vec3& vv):v1(vd.v),v2(vv),d(vd.d){}

  //operator Vec3();
};

inline VDMulVadd operator+(const VDMul vd,const Vec3& vv)
{
  return VDMulVadd(vd,vv);
}

#endif
