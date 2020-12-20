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

#include "Edge.h"

using std::make_pair;

/*!
  Construct edge from 2 points. Set triangle pointers to NULL
  
  \param id1 the id of p1
  \param id2 the id of p2
  \param p1
  \param p2
*/
Edge::Edge(int id1,int id2,const Vec3& p1,const Vec3& p2) : AEdge(p1,p2)
{
  m_id1=id1;
  m_id2=id2;
  m_t1=NULL;
  m_t2=NULL;
}

void Edge::applyForce(const Vec3 &f)
{
  if ((m_t1 != NULL) && (m_t2 != NULL))
  {
    const Vec3 halfF = f*(0.5);
    m_t1->applyForce(halfF);
    m_t2->applyForce(halfF);
  }
  else if (m_t1 != NULL)
  {
    m_t1->applyForce(f);
  }
  else if (m_t2 != NULL)
  {
    m_t2->applyForce(f);
  }
}

/*!
  Construct edge from 2 points and 1 triangle pointer. Set other triangle pointer to NULL
  
  \param id1 the id of p1
  \param id2 the id of p2
  \param p1
  \param p2
  \param t1
*/
Edge::Edge(int id1,int id2,const Vec3& p1,const Vec3& p2,Triangle* t1) : AEdge(p1,p2)
{
  m_id1=id1;
  m_id2=id2;
  m_t1=t1;
  m_t2=NULL;
}

/*!
  Construct edge from 2 points and 2 triangle pointers.
  
  \param id1 the id of p1
  \param id2 the id of p2
  \param p1
  \param p2
  \param t1
  \param t2
*/
Edge::Edge(int id1,int id2,const Vec3& p1,const Vec3& p2,Triangle* t1,Triangle* t2) : AEdge(p1,p2)
{
  m_id1=id1;
  m_id2=id2;
  m_t1=t1;
  m_t2=t2;
}

/*!
  Check if any of the adjacent triangles (if there are any) has contact, i.e
  the perpendicular line from the supporting plane to the point hits the triangle,
  an thus makes the edge contact invalid

  \param P the point
*/
bool Edge::isValidContact(const Vec3& P) const
{
  bool t1,t2,is_valid;

  // check 1st triangle 
  if(m_t1!=NULL) { 
    t1=(m_t1->dist(P)).first;
  } else {
    t1=false;
  }
  // check 2nd triangle
  if(m_t2!=NULL) { 
    t2=(m_t2->dist(P)).first;
  } else {
    t2=false;
  }
  is_valid=(!t1) && (!t2);
  return is_valid;
}

/*!
  Get min. corner of axis-aligned bounding box.
*/
Vec3 Edge::getBoundingBoxMin() const
{
  return cmin(m_p0,m_p1);
}
 
/*!
  Get min. corner of axis-aligned bounding box.
*/
Vec3 Edge::getBoundingBoxMax() const
{
  return cmax(m_p0,m_p1);
}

/*!
  get unit direction vector between a point and the closest
  point along the supporting line of the edge (pointing away 
  from the edge). 

  \param p the point
  \warning does not check if the closest point is actually within the edge or 
  that the point is not on the line (potential div by 0) 
*/
Vec3 Edge::getDirectionFromPoint(const Vec3& p) const
{
  Vec3 dir;

  Vec3 v=m_p1-m_p0;
  Vec3 vu=v.unit();
  double d=((p-m_p0)*vu);
  dir=(p-m_p0)-d*vu;

  return dir.unit();
}

/*!
  Move one of the corners. The identifier for the corner is the global node id.
  If the edge doesn't contain the node with the id , do nothing.

  \param id the global id of the node to be moved
  \param d the amount of movement
*/
void Edge::moveNode(int id,const Vec3& d)
{
  if (id==m_id1){
    m_p0+=d;
  } else if (id==m_id2){
    m_p1+=d;
  } else {
    std::cerr << "trying to move node not in edge!" << std::endl;
  }
}

/*!
  Translate whole edge

  \param d the amount of movement
*/
void Edge::move(const Vec3& d)
{
  m_p0+=d;
  m_p1+=d;
}

void Edge::rotate(const Vec3& origin ,const Vec3& axis,double angle)
{
    m_p0.rotateBy(origin,axis,angle);
    m_p1.rotateBy(origin,axis,angle);
}

/*!
  output Edge to ostream
*/
ostream& operator << (ostream& ost,const Edge& E)
{
  ost << "Edge: (" << E.m_p0 << ") - (" << E.m_p1 << ")";
  if(E.m_t1!=NULL) ost << " in : " << *(E.m_t1);
  if(E.m_t2!=NULL) ost << " in : " << *(E.m_t2);
  
  return ost;
}
