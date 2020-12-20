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

#ifndef __ABCDAMPING_IGP_H
#define __ABCDAMPING_IGP_H

// -- project includes --
#include "Model/DampingIGP.h"
#include "Foundation/vec3.h"

/*!
  Interaction group parameters for ABCDamping
*/
class ABCDampingIGP : public CDampingIGP 
{
 protected:
  Vec3 m_pos,m_normal;
  double m_c1;

 public:
  ABCDampingIGP(){};
  ABCDampingIGP(const string&,const string&,double,double,int,const Vec3&,const Vec3&,const Vec3&,double);

  virtual void  packInto(CVarMPIBuffer*) const;

  void setPos(const Vec3& p){m_pos=p;};
  Vec3 getPos(){return m_pos;};
  void setNormal(const Vec3& n){m_normal=n;};
  Vec3 getNormal(){return m_normal;};
  void setC1(double d){m_c1=d;};
  double getC1(){return m_c1;};
};

ABCDampingIGP* extractABCDampingIGP(AMPIBuffer*);

#endif // __ABCDAMPING_IGP_H
