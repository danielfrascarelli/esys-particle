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

#include "vec3.h"

//the error...
VecErr::VecErr(const string& m):MError(m)
{
  message.insert(0,"Vec3 "); 
}

// constructors
Vec3::Vec3()
{
  data[0]=0;
  data[1]=0;
  data[2]=0;
}

Vec3::Vec3(double a,double b,double c)
{
  data[0]=a;
  data[1]=b;
  data[2]=c;
}

Vec3::Vec3(const Vec3& rhs)
{
  data[0]=rhs.data[0];
  data[1]=rhs.data[1];
  data[2]=rhs.data[2];
}

// operators

Vec3& Vec3::operator=(const Vec3& rhs)
{
  data[0]=rhs.data[0];
  data[1]=rhs.data[1];
  data[2]=rhs.data[2];
  return *this;
}

Vec3& Vec3::operator-=(const Vec3& rhs)
{
  data[0]-=rhs.data[0];
  data[1]-=rhs.data[1];
  data[2]-=rhs.data[2];
  return *this;
}

Vec3& Vec3::operator+=(const Vec3& rhs)
{
  data[0]+=rhs.data[0];
  data[1]+=rhs.data[1];
  data[2]+=rhs.data[2];
  return *this;
}

Vec3 Vec3::operator+(const Vec3& rhs) const
{
  return Vec3(data[0]+rhs.data[0], data[1]+rhs.data[1], data[2]+rhs.data[2]);
}

Vec3 Vec3::operator-(const Vec3& rhs) const
{
  return Vec3(data[0]-rhs.data[0], data[1]-rhs.data[1], data[2]-rhs.data[2]);
}

double Vec3::operator*(const Vec3& rhs) const
{
  return data[0]*rhs.data[0]+data[1]*rhs.data[1]+data[2]*rhs.data[2];
}

Vec3 Vec3::operator*(double s) const 
{ 
   return Vec3(data[0]*s,data[1]*s,data[2]*s) ; 
} 

Vec3 Vec3::operator/(double s) const 
{ 
   return Vec3(data[0]/s,data[1]/s,data[2]/s) ; 
} 


// vector product
// 9 Flops ( 6 mult, 3 sub ) 
Vec3 cross(const Vec3& lhs,const Vec3& rhs)
{
  return Vec3(lhs.data[1]*rhs.data[2]-lhs.data[2]*rhs.data[1],
	      lhs.data[2]*rhs.data[0]-lhs.data[0]*rhs.data[2],
	      lhs.data[0]*rhs.data[1]-lhs.data[1]*rhs.data[0]);
}
  
Vec3 operator*(double f,const Vec3& rhs)
{
  return Vec3(f*rhs.data[0], f*rhs.data[1], f*rhs.data[2]);
}


// euclidian norm
// 6 Flops ( 3 mult, 2 add, 1 sqrt )
double Vec3::norm() const
{
  return sqrt(data[0]*data[0]+data[1]*data[1]+data[2]*data[2]);
}

// square of the euclidian norm
// 5 Flops ( 3 mult, 2 add)
double Vec3::norm2() const
{
  return data[0]*data[0]+data[1]*data[1]+data[2]*data[2];
}

// returns unit vector in direction of the original vector
// 9 Flops ( 3 mult, 2 add, 3 div, 1 sqrt ) 
Vec3 Vec3::unit() const
{
  Vec3 res(data[0],data[1],data[2]);
  return res/norm();
}

// per element min/max
Vec3 cmax(const Vec3& v1,const Vec3& v2)
{
  Vec3 res;
  res.data[0]=v1.data[0]>v2.data[0] ? v1.data[0] : v2.data[0];
  res.data[1]=v1.data[1]>v2.data[1] ? v1.data[1] : v2.data[1];
  res.data[2]=v1.data[2]>v2.data[2] ? v1.data[2] : v2.data[2];
  return res;
}

Vec3 cmin(const Vec3& v1,const Vec3& v2)
{
  Vec3 res;
  res.data[0]=v1.data[0]<v2.data[0] ? v1.data[0] : v2.data[0];
  res.data[1]=v1.data[1]<v2.data[1] ? v1.data[1] : v2.data[1];
  res.data[2]=v1.data[2]<v2.data[2] ? v1.data[2] : v2.data[2];
  return res;
}

// save version, throws exception if norm()==0
Vec3 Vec3::unit_s() const
{
  double n=norm();
  if(n==0) throw VecErr("norm() of data[2]ero-vector"); 
  Vec3 res(data[0],data[1],data[2]);
  return res/n;
}



double Vec3::max() const
{
  double m;

  m=( data[0]>data[1] ? data[0] : data[1] );
  m=( m>data[2] ? m : data[2] );

  return m;
}

double Vec3::min() const
{
  double m;

  m=( data[0]<data[1] ? data[0] : data[1] );
  m=( m<data[2] ? m : data[2] );

  return m;
}
  
bool Vec3::operator==(const Vec3& V)
{
  return((data[0]==V.data[0])&&(data[1]==V.data[1])&&(data[2]==V.data[2]));
}

bool Vec3::operator!=(const Vec3& V)
{
  return((data[0]!=V.data[0])||(data[1]!=V.data[1])||(data[2]!=V.data[2]));
}

// in/output

ostream& operator << (ostream& ostr,const Vec3& V)
{
  ostr << "( " ;
  ostr << V.data[0] << " , " << V.data[1] << " , " << V.data[2] <<" )";
  return ostr;
}

istream& operator >> (istream& istr,Vec3& V)
{
  istr >> V.data[0];
  istr >> V.data[1];
  istr >> V.data[2];
  
  return istr;
}

// n+1-ary operators

Vec3::Vec3(const VDMulVadd& v)
{
  mul_add_and_assign(&v.v1,&v.v2,v.d);
}

Vec3& Vec3::operator=(const VDMulVadd& v)
{
  mul_add_and_assign(&v.v1,&v.v2,v.d);
  return *this;
}
  
Vec3::Vec3(const VDMul& v)
{
  mul_and_assign(&v.v,v.d);
}

Vec3& Vec3::operator=(const VDMul& v)
{
  mul_and_assign(&v.v,v.d);
  return *this;
}
  
void Vec3::mul_add_and_assign(const Vec3* v1,const Vec3* v2,const double& d)
{
  data[0]=v1->data[0]*d+v2->data[0];
  data[1]=v1->data[1]*d+v2->data[1];
  data[2]=v1->data[2]*d+v2->data[2];
}

void Vec3::mul_and_assign(const Vec3* v,const double& d)
{
  data[0]=v->data[0]*d;
  data[1]=v->data[1]*d;
  data[2]=v->data[2]*d;
}
