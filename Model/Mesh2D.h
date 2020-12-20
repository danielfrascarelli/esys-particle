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

#ifndef __MESH2D_H
#define __MESH2D_H

// -- Project includes --
#include "Geometry/Edge2D.h"
#include "Geometry/Corner2D.h"
#include "Model/MeshData2D.h"
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

//--- TML includes ---
#include "tml/comm/comm.h"

/*!
  A class for a 2D "mesh", i.e. a collection of lines
  connected by corners in a 2D plane which can interact with 
  particles in 2D. The Main purpose is for coupling with 2D FEM
  (Finley) simulations via escript.
*/
class Mesh2D
{
 private:
  vector<Edge2D> m_edges;
  vector<Corner2D> m_corners;
  map<int,int> m_corner_by_id;
  multimap<int,Edge2D*> m_edge_by_node_id;
  map<int,int> m_edge_index_by_id;

 public:
  // types 
  typedef vector<Edge2D>::iterator edge_iterator;
  typedef vector<Corner2D>::iterator corner_iterator;

  // functions
  Mesh2D();
  virtual ~Mesh2D(){};
  void LoadMesh(const vector<MeshNodeData2D>&,const vector<MeshEdgeData2D>&);
  void moveNode(int,const Vec3&);
  void translateBy(const Vec3 &translation);

  edge_iterator edges_begin(){return m_edges.begin();};
  edge_iterator edges_end(){return m_edges.end();};
  corner_iterator corners_begin(){return m_corners.begin();};
  corner_iterator corners_end(){return m_corners.end();};
  Edge2D* getEdgeById(int);
  Corner2D* getCornerById(int);

  void zeroForces();
  virtual void writeCheckPoint(ostream&,const string&) const;
  virtual void loadCheckPoint(istream&);
  
  // edge data access functions
  template <typename P> void forAllEdgesGet(P&,typename P::value_type (Edge2D::*rdf)() const);
  template <typename P> vector<pair<int,P> > forAllEdgesGetIndexed(P (Edge2D::*rdf)() const);
};

#include "Mesh2D.hpp"

#endif // __MESH2D_H
