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

#ifndef __CORNER_H
#define __CORNER_H

//-- Project includes --
#include "Foundation/vec3.h"
#include "Geometry/Triangle.h"
#include "Geometry/Edge.h"

//-- STL includes --
#include <vector>

using std::vector;

/*!
  \class Corner
  \brief Class representing the corner of a polygon

  \author Steffen Abe
  $Revision$
  $Date$
*/
class Corner
{
 private:
  Vec3 m_p;
  Vec3 m_old_pos;
  vector<Edge*> m_edges;
  vector<Triangle*> m_triangles;
  int m_id;
  int m_tag;

 public:
  Corner(const Vec3&,int,int);

  void addEdge(Edge*);
  void addTriangle(Triangle*);
  double sep(const Vec3&) const;
  //  pair<bool,double> dist(const Vec3&) const ; // signed separation according to direction of the normal
  bool isValidContact(const Vec3&) const;
  Vec3 getDirectionFromPoint(const Vec3&) const;
  void move(const Vec3&);
  void rotate(const Vec3&,const Vec3&,double);
  Vec3 getPos()const {return m_p;};
  void setPos(const Vec3 &p) {m_p = p;}
  void applyForce(const Vec3 &f);
  int getID() const {return m_id;};
  int getTag() const {return m_tag;};
  
  double getDistMoved() {return (m_old_pos-m_p).norm();};
  void resetOldPos(){m_old_pos=m_p;};
};

#endif // __CORNER_H
