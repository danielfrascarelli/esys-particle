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

#ifndef __EDGE2D_H
#define __EDGE2D_H

//-- Project includes --
#include "Foundation/vec3.h"
#include "Foundation/Matrix3.h"
#include "Geometry/AEdge.h"

//-- STL includes --
#include <utility>
using std::pair;
using std::make_pair;

//-- IO includes --
#include <iostream>
using std::ostream;

/*!
  \class Edge2D
  \brief class for edge in 2D "mesh"

  \author Steffen Abe
  $Revision$
  $Date$
*/
class Edge2D : public AEdge
{
 public: // types
  typedef Vec3 (Edge2D::* VectorFieldFunction)() const;
  typedef double (Edge2D::* ScalarFieldFunction)() const;
  
 private:
  Vec3 m_normal;
  Vec3 m_force;
  int m_id0,m_id1; // corner ids
  int m_edge_id,m_tag;

 public:
  Edge2D(int,int,const Vec3&,const Vec3&,int,int);
  void moveNode(int,const Vec3&);

  inline int getID(){return m_edge_id;};
  inline void applyForce(const Vec3& f){m_force+=f;};
  inline void zeroForce(){m_force=Vec3(0.0,0.0,0.0);};
  Vec3 getNormal() const {return m_normal;};
  Vec3 toGlobal(const Vec3&);
  Vec3 toLocal(const Vec3&);

  // get id/pos pairs for each node -> mainly for checkpointing
  pair<int,Vec3> getP0()const{return make_pair(m_id0,m_p0);};
  pair<int,Vec3> getP1()const{return make_pair(m_id1,m_p1);};
  
  // access functions
  static VectorFieldFunction getVectorFieldFunction(const string&);
  static ScalarFieldFunction getScalarFieldFunction(const string&);

  Vec3 getForce() const {return m_force;};
  Vec3 getForceDensity() const {return m_force/((m_p1-m_p0).norm());};
  double getPressure() const;

  //! output for debugging purposes
  friend ostream& operator<<(ostream&,const Edge2D&); 
  void print();
};


#endif // __EDGE2D_H
