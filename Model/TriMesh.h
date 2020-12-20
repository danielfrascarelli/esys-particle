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

#ifndef __TRIMESH_H
#define __TRIMESH_H

//--- TML includes ---
#include "tml/comm/comm.h"

// -- Project includes --
#include "Geometry/Triangle.h"
#include "Geometry/Edge.h"
#include "Geometry/Corner.h"
#include "Model/MeshData.h"
#include "Foundation/vec3.h"

// -- STL includes --
#include <vector>
#include <map>
#include <string>

using std::vector;
using std::multimap;
using std::map;
using std::string;

// -- IO includes ---
#include <iostream>

using std::ostream;


/*!
  \class TriMesh
  \brief class for a triangle mesh

  \author Steffen Abe
  $Revision$
  $Date$
*/
class TriMesh
{
 private:
  vector<Triangle> m_triangles;
  vector<Edge> m_edges;
  vector<Corner> m_corners;
  multimap<int,Triangle*> m_triangle_by_node_id;
  multimap<int,Edge*> m_edge_by_node_id;
  map<int,int> m_corner_by_id;
  
  map<int,int> m_tri_index_by_id;

 public:
  // types 
  typedef vector<Triangle>::iterator triangle_iterator;
  typedef vector<Edge>::iterator edge_iterator;
  typedef vector<Corner>::iterator corner_iterator;
 
  // functions
  TriMesh();

  virtual ~TriMesh()
  {
  }
  
  void LoadMesh(const vector<MeshNodeData>&,const vector<MeshTriData>&);
  void moveNode(int,const Vec3&);
  void translateBy(const Vec3 &translation);
  void rotateBy(const Vec3&,const Vec3 &axis, double angle);
  triangle_iterator triangles_begin(){return m_triangles.begin();};
  triangle_iterator triangles_end(){return m_triangles.end();};
  edge_iterator edges_begin(){return m_edges.begin();};
  edge_iterator edges_end(){return m_edges.end();};
  corner_iterator corners_begin(){return m_corners.begin();};
  corner_iterator corners_end(){return m_corners.end();};
  Triangle* getTriangleById(int);
  bool hasMovedBy(double);
  void resetCurrentDisplacement();

  void zeroForces();
  virtual void writeCheckPoint(ostream&,const string&) const;
  virtual void loadCheckPoint(istream&);

  // triangle data access functions
  template <typename P> void forAllTrianglesGet(P&,typename P::value_type (Triangle::*rdf)() const);
  template <typename P> vector<pair<int,P> > forAllTrianglesGetIndexed(P (Triangle::*rdf)() const);
};

#include "TriMesh.hpp"

#endif // __TRIMESH_H
