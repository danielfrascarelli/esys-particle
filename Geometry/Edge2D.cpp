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


#include "Foundation/console.h"
#include "Geometry/Edge2D.h"

// --- System includes
#include <cmath>

using std::fabs;
using std::endl;
using std::make_pair;

/*!
  construct Edge from corner coordinates. Make sure Z=0 for all corners

  \param id0 id of the first corner
  \param id1 id of the 2nd corner
  \param v0 first corner
  \param v1 second corner
  \param ed_id edge id
  \param tag edge tag

*/  
Edge2D::Edge2D(int id0,int id1,const Vec3& v0,const Vec3& v1,int ed_id,int tag) : AEdge(v0,v1)
{
  m_id0=id0;
  m_id1=id1;
  m_p0.Z()=0.0;
  m_p1.Z()=0.0;
  m_edge_id=ed_id;
  m_tag=tag;
  m_force=Vec3(0.0,0.0,0.0);
  m_normal=cross(Vec3(0.0,0.0,1.0),m_p1-m_p0).unit();
}

/*!
  Move one of the corners. The identifier for the corner is the global node id.
  If the node with the id is not in the edge do nothing.

  \param id the global id of the node to be moved
  \param d the amount of movement
*/
void Edge2D::moveNode(int id,const Vec3& d)
{
  if(id==m_id0){
    m_p0+=d;
    m_normal=cross(Vec3(0.0,0.0,1.0),m_p1-m_p0).unit();
  } else if (id==m_id1) {
    m_p1+=d;
    m_normal=cross(Vec3(0.0,0.0,1.0),m_p1-m_p0).unit();
  } 
}

/*!
  convert point in local coords to global coords

  \param p the point in local coords
*/
Vec3 Edge2D::toGlobal(const Vec3& p)
{
  Vec3 res;

  res=m_p0+(m_p1-m_p0)*p.X()+m_normal*p.Y();

  return res;
}

/*!
  convert point in global coords into local (x,y,0) coords

  \param p the point in global coords 
*/
Vec3 Edge2D::toLocal(const Vec3& p)
{
  double x;
  double y;

  x=((p-m_p0)*(m_p1-m_p0))/(m_p1-m_p0).norm2();
  y=(p-m_p0)*m_normal;

  return Vec3(x,y,0.0);
}

/*!
  Get the Edge2D member function which returns a vector field of a given name. Returns
  NULL if a field with that name doesn't exist.

  \param name the name of the field 
*/
Edge2D::VectorFieldFunction Edge2D::getVectorFieldFunction(const string& name)
{
  Edge2D::VectorFieldFunction f;

  if(name=="force"){
    f=&Edge2D::getForce;
  } if(name=="forcedensity"){
    f=&Edge2D::getForceDensity;
  } else {
    f=NULL;
    cerr << "ERROR - invalid name for edge vector access function" << endl; 
  }

  return f;
}

/*!
  Get the Edge2D member function which returns a scalar field of a given name. Returns
  NULL if a field with that name doesn't exist.

  \param name the name of the field 
*/
Edge2D::ScalarFieldFunction Edge2D::getScalarFieldFunction(const string& name)
{
  Edge2D::ScalarFieldFunction f;

  if(name=="pressure"){
    f=&Edge2D::getPressure;
  } else {
    f=NULL;
    cerr << "ERROR - invalid name for edge scalar access function" << endl; 
  }

  return f;
}


/*!
  Get pressure on the edge from interaction forces
*/
double Edge2D::getPressure() const
{
  // calculate normal force
  double F_n=m_force*m_normal;
  // retrun force per area (length)
  return F_n/(m_p1-m_p0).norm();    
}

/*!
  output Edge2D to ostream
*/
ostream& operator<<(ostream& ost,const Edge2D& T)
{
  ost << "Edge2D: (" << T.m_p0 << ") - (" << T.m_p1 << ") Normal: (" << T.m_normal << ")"; 
  return ost;
}

/*!
  output Edge2D to cout
*/
void Edge2D::print()
{
  cout << "Edge2D: (" << m_p0 << ") - (" << m_p1 << ") Normal: (" << m_normal << ")"; 
}
