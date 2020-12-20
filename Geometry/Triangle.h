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

#ifndef __TRIANGLE_H
#define __TRIANGLE_H


//-- Project includes --
#include "Foundation/vec3.h"
#include "Foundation/Matrix3.h"

//-- STL includes --
#include <utility>
using std::pair;
using std::make_pair;

//-- IO includes --
#include <iostream>
using std::ostream;

//! exception class for Triangle
class TriangleError
{
public:
  TriangleError(){};
};


/*!
  \class Triangle
  \brief Class representing a Triangle

  \author Steffen Abe
  $Revision$
  $$Date$

*/
class Triangle
{
 public: // types
  typedef Vec3 (Triangle::* VectorFieldFunction)() const;
  typedef double (Triangle::* ScalarFieldFunction)() const;

 private:
  Matrix3 m_invtrans;
  Matrix3 m_trans;
  Vec3 m_p0,m_p1,m_p2;
  Vec3 m_normal;
  Vec3 m_force;
  int m_id0,m_id1,m_id2;
  int m_tri_id,m_tag;

  double EdgeSep(const Vec3&, const Vec3& ,const Vec3& ) const;
  
 public:
  Triangle(int,int,int,const Vec3&,const Vec3&,const Vec3&,int,int);

  double sep(const Vec3&) const;
  pair<bool,double> dist(const Vec3&) const ; // signed separation according to direction of the normal
  Vec3 getBoundingBoxMin() const; 
  Vec3 getBoundingBoxMax() const; 
  Vec3 getNormal() const {return m_normal;};
  Vec3 toGlobal(const Vec3&);
  Vec3 toLocal(const Vec3&);
  bool containsEdge(const Vec3&,const Vec3&) const;
  void moveNode(int,const Vec3&);
  void move(const Vec3&);
  void rotate(const Vec3&,const Vec3&,double);
  inline int getID() const {return m_tri_id;};
  inline int getTag() const  {return m_tag;}; 
  inline void applyForce(const Vec3& f){m_force+=f;};
  inline void zeroForce(){m_force=Vec3(0.0,0.0,0.0);};

  // get id/pos pairs for each node -> mainly for checkpointing
  pair<int,Vec3> getP0()const{return make_pair(m_id0,m_p0);};
  pair<int,Vec3> getP1()const{return make_pair(m_id1,m_p0+m_p1);};
  pair<int,Vec3> getP2()const{return make_pair(m_id2,m_p0+m_p2);};

  inline int getPid0() const{return m_id0;};
  inline int getPid1() const{return m_id1;};
  inline int getPid2() const{return m_id2;};

  // access functions
  static VectorFieldFunction getVectorFieldFunction(const string&);
  static ScalarFieldFunction getScalarFieldFunction(const string&);

  Vec3 getForce()const {return m_force;};
  double getPressure() const;

  //! output for debugging purposes
  friend ostream& operator<<(ostream&,const Triangle&); 
};

#endif //__TRIANGLE_H
