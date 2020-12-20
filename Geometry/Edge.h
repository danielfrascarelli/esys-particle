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

#ifndef __EDGE_H
#define __EDGE_H

//-- Project includes --
#include "Foundation/vec3.h"
#include "Geometry/Triangle.h"
#include "Geometry/AEdge.h"

//-- STL includes --
#include <utility>

using std::pair;
using std::make_pair;

/*!
  \class Edge
  \brief Class representing the edge of a polygon

  \author Steffen Abe
  $Revision$
  $Date$
*/
class Edge : public AEdge
{
 private:
  Triangle *m_t1,*m_t2;
  int m_id1,m_id2;

 public:
  Edge(int,int,const Vec3&,const Vec3&);
  Edge(int,int,const Vec3&,const Vec3&,Triangle*);
  Edge(int,int,const Vec3&,const Vec3&,Triangle*,Triangle*);

  bool isValidContact(const Vec3&) const;
  Vec3 getBoundingBoxMin() const; 
  Vec3 getBoundingBoxMax() const; 
  Vec3 getDirectionFromPoint(const Vec3&) const;
  void moveNode(int,const Vec3&);
  void move(const Vec3&);
  void rotate(const Vec3&,const Vec3&,double);
  void applyForce(const Vec3 &f);

  pair<int,int> getIDs() const {return make_pair(m_id1,m_id2);};

  friend ostream& operator << (ostream&,const Edge&);
};
#endif // __EDGE_H
