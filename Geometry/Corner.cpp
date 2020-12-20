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

#include "Corner.h"

/*!
  constructor

  \param pos the position of the corner
  \param id the node id
  \param tag the node tag
*/
Corner::Corner(const Vec3& pos,int id, int tag)
{
  m_p=pos;
  m_old_pos=m_p;
  m_id=id;
  m_tag=tag;
}
 
/*!
  add Edge to corner
  
  \param edge a pointer to the edge
*/ 
void Corner::addEdge(Edge* edge)
{
  m_edges.push_back(edge);
}

/*!
  add Triangle to Corner

  \param triangle a pointer to the triangle
*/
void Corner::addTriangle(Triangle* triangle)
{
  m_triangles.push_back(triangle);
}

void Corner::applyForce(const Vec3 &f)
{
  if (m_triangles.size() > 0)
  {
    const Vec3 portionF = f*(1.0/m_triangles.size());
    for (
      std::vector<Triangle *>::iterator it = m_triangles.begin();
      it != m_triangles.end();
      ++it
    ){
      (*it)->applyForce(portionF);
    }
  }
}

/*!
  get distance between corner and point

  \param p the point
*/
double Corner::sep(const Vec3& p) const
{
  return (m_p-p).norm();
}

/*!
  check if the contact between a particle at a point and the 
  corner is valid or if there is a contact between the particle
  and any of the adjacent edges or triangles

  \param p the center of the particle 
*/
bool Corner::isValidContact(const Vec3& p) const
{ 
  bool res=true;
  
  // check vs. Triangles
  vector<Triangle*>::const_iterator titer=m_triangles.begin();
  while((titer!=m_triangles.end())&&(res)){
    res=!((*titer)->dist(p).first);
    titer++;
  }

  // if no contact with triangle, check vs. Edges
  if(res){
    vector<Edge*>::const_iterator eiter=m_edges.begin();
    while((eiter!=m_edges.end())&&(res)){
      res=!((*eiter)->dist(p).first);
      eiter++;
    }
  }
  return res;
}

/*!
  get the unit direction vector between a point and the corner

  \param p the point
*/
Vec3 Corner::getDirectionFromPoint(const Vec3& p) const
{
  return (p-m_p).unit();
}

/*!
  move the corner

  \param d the displacement
*/
void Corner::move(const Vec3& d)
{
  m_p+=d;
}


void Corner::rotate(const Vec3& origin,const Vec3& axis,double angle)
{
    m_p.rotateBy(origin,axis,angle);
}