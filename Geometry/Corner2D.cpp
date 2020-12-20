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

#include "Corner2D.h"
#include "console.h"

/*!
  constructor, make sure Z=0;

  \param pos the position of the corner
  \param id the node id
*/
Corner2D::Corner2D(const Vec3& pos, int id)
{
  m_p=pos;
  m_p.Z()=0.0;
  m_id=id;
}


/*!
  add Edge to corner
  
  \param edge a pointer to the edge
*/ 
void Corner2D::addEdge(Edge2D* edge)
{
  m_edges.push_back(edge);
}

/*!
  get distance between corner and point

  \param p the point
*/
double Corner2D::sep(const Vec3& p) const
{
  return (m_p-p).norm();
}


/*!
  check if the contact between a particle at a point and the 
  corner is valid or if there is a contact between the particle
  and any of the adjacent edges or triangles

  \param p the center of the particle 
*/
bool Corner2D::isValidContact(const Vec3& p) const
{
  bool res=true;

  // check vs. Edges
  vector<Edge2D*>::const_iterator eiter=m_edges.begin();
  while((eiter!=m_edges.end())&&(res)){
    res=!((*eiter)->dist(p).first);
    eiter++;
  }

  return res;
}

/*!
  get the unit direction vector between a point and the corner

  \param p the point
*/
Vec3 Corner2D::getDirectionFromPoint(const Vec3& p) const
{
  return (p-m_p).unit();
}

/*!
  move the corner, make sure Z=0;

  \param d the displacement
*/
void Corner2D::move(const Vec3& d)
{
  m_p+=d;
  m_p.Z()=0.0;
}

/*!
  get normal of an edge

  \param idx which egde (1,2) to get
*/
Vec3 Corner2D::getEdgeNormal(int idx) const
{
  Vec3 res;

  if((idx==1) && (m_edges.size()>=1)) {
    res=(m_edges[0])->getNormal();
  } else if ((idx==2) && (m_edges.size()>=2)) {
    res=(m_edges[1])->getNormal();
  } else {
    console.Error() << "Error in Corner2D::getEdgeNormal: idx=" << idx << " nr. of edges: " << m_edges.size() << "\n";
  }

  return res;
}

/*!
  apply force to one of the attached edges

  \param idx which egde (1,2)
  \param F the force
*/
void Corner2D::applyForceToEdge(int idx,const Vec3& f)
{
  if((idx==1) && (m_edges.size()>=1)) {
    (m_edges[0])->applyForce(f);
  } else if ((idx==2) && (m_edges.size()>=2)) {
    (m_edges[1])->applyForce(f);
  } else {
    console.Error() << "Error in Corner2D::applyForceToEdge : idx=" << idx << " nr. of edges: " << m_edges.size() << "\n";
  }
}
