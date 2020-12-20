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

#ifndef __CORNER2D_H
#define __CORNER2D_H

//-- Project includes --
#include "Foundation/vec3.h"
#include "Geometry/Edge2D.h"

//-- STL includes --
#include <vector>

using std::vector;


/*!
  \class Corner2D 
  \brief Class representing the corner in a 2D "mesh"

  \author Steffen Abe
  $Revision$
  $Date$
*/
class Corner2D
{
 private:
  Vec3 m_p;
  vector<Edge2D*> m_edges;
  int m_id;

 public:
  Corner2D(const Vec3&,int);

  void addEdge(Edge2D*);
  double sep(const Vec3&) const;
  //  pair<bool,double> dist(const Vec3&) const ; // signed separation according to direction of the normal
  bool isValidContact(const Vec3&) const;
  Vec3 getDirectionFromPoint(const Vec3&) const;
  void move(const Vec3&);
  Vec3 getPos()const {return m_p;};
  int getID() const {return m_id;};
  int getNEdges() const {return m_edges.size();};
  Vec3 getEdgeNormal(int) const;
  void applyForceToEdge(int,const Vec3&);
};

#endif // __CORNER2D_H
