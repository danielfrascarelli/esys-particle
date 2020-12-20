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

#include "Triangle.h"
#include "Foundation/console.h"

// --- System includes
#include <cmath>

using std::fabs;
using std::endl;
using std::make_pair;

/*!
  Construct triangle from the corner coordinates. It is assumed that the corners
  are given anticlockwise and the normal is calculated accordingly.

  \param id0 id of the first corner
  \param id1 id of the 2nd corner
  \param id2 id of the 3rd corner
  \param v0 first corner
  \param v1 second corner
  \param v2 third corner
  \param tri_id triangle id
  \param tag triangle tag
*/
Triangle::Triangle(int id0,int id1,int id2,const Vec3& v0,const Vec3& v1,const Vec3& v2,int tri_id,int tag)
{
  m_id0=id0;
  m_id1=id1;
  m_id2=id2;
  m_p0=v0; 
  m_p1=v1-m_p0;
  m_p2=v2-m_p0;
  m_tri_id=tri_id;
  m_tag=tag;
  m_force=Vec3(0.0,0.0,0.0);
  try{
    m_normal=cross(m_p2,m_p1).unit_s();
    m_trans=Matrix3(m_p1,m_p2,m_normal);
    m_invtrans=m_trans.inv();
  }
  catch(VecErr V){ // if unit_s thows exception -> v1||v2
    throw TriangleError();
  }
}

/*!
  calculate distance between an edge and a point

  \param p0 point 1 of the edge
  \param p1 point 2 of the edge
  \param p the point
*/
double Triangle::EdgeSep(const Vec3& p0,const Vec3& p1,const Vec3& p) const
{
  double sep;

  Vec3 v=p1-p0;
  Vec3 vu=v.unit();
  double d=((p-p0)*vu);
  if((d>0.0)&(d<v.norm())){
    sep=((p-p0)-d*vu).norm();
  }else{
    sep=-1;
  }
  return sep;
}

/*!
  Get distance between point and the triangle. 

  \param p the point
*/
double Triangle::sep(const Vec3& p) const
{
  double sep=-1;

  // transform point to triangle local coord system
  Vec3 p_local=(p-m_p0)*m_invtrans;
  // check if closest point is in triangle
  if((p_local.X()>=0.0) && (p_local.Y()>=0.0) &&(p_local.X()+p_local.Y()<=1.0)){
    sep=fabs((p-m_p0)*m_normal);
  } else { // need to check distance to edges/corners 
    double d1=EdgeSep(m_p0,m_p0+m_p1,p);
    double d2=EdgeSep(m_p0,m_p0+m_p2,p);
    double d3=EdgeSep(m_p0+m_p1,m_p0+m_p2,p);
    // find the minimum separation != -1 (messy)
    if(d1>0.0){
      if(d2>0.0){
	sep=(d1<d2) ? d1 : d2;
	if(d3>0.0){
	  sep=(d3<sep) ? d3 : sep; 
	}
      } else if(d3>0){
	sep=(d1<d3) ? d1 : d3;
      } else {
	sep=d1;
      }
    } else if (d2>0){
      if (d3>0){
	sep=(d2<d3) ? d2 : d3;
      } else {
	sep=d2;
      }
    } else {
      sep=d3;
    }
    if(sep==-1.0){ // no edge-> get corner dist
      d1=(p-m_p0).norm();
      d2=(p-(m_p0+m_p1)).norm();
      d3=(p-(m_p0+m_p2)).norm();
      sep=(d1<d2) ? d1 : d2;
      sep=(sep<d3) ? sep : d3;
    }
  }

  return sep;
}

/*!
  Get the signed distance between a point and the triangle.
  If the closest point on the supporting plane is outside the 
  triangle, the first component of the return value is "false", 
  otherwise "true"


  \param p the point
*/
pair<bool,double> Triangle::dist(const Vec3& p) const 
{
  bool is_in=false;
  double dist=0.0;

  // transform point to triangle local coord system
  Vec3 p_local=(p-m_p0)*m_invtrans;
  //cout << "p_local : " << p_local << endl;
  // check if closest point is in triangle
  if((p_local.X()>=0.0) && (p_local.Y()>=0.0) &&(p_local.X()+p_local.Y()<=1.0)){
    dist=(p-m_p0)*m_normal;
    is_in=true;
  } else { 
    is_in=false;
  }
  return make_pair(is_in,dist);
}

/*!
  Get min. corner of axis-aligned bounding box of the triangle.
*/
Vec3 Triangle::getBoundingBoxMin() const
{
  Vec3 v_min=cmin(m_p0,m_p0+m_p1);
  v_min=cmin(v_min,m_p0+m_p2);

  return v_min;
}

/*!
  Get max. corner of axis-aligned bounding box of the triangle.
*/
Vec3 Triangle::getBoundingBoxMax() const
{  
  Vec3 v_max=cmax(m_p0,m_p0+m_p1);
  v_max=cmax(v_max,m_p0+m_p2);

  return v_max;
}

/*!
  check if an edge given by 2 points is in the triangle

  \param p1 
  \param p2
*/
bool Triangle::containsEdge(const Vec3& p1,const Vec3& p2) const
{
  bool p1_in=((p1==m_p0)||(p1==m_p0+m_p1) || (p1==m_p0+m_p2));
  bool p2_in=((p2==m_p0)||(p2==m_p0+m_p1) || (p2==m_p0+m_p2));
  return ((p1!=p2) && p1_in && p2_in);
}

/*!
  Move one of the corners. The identifier for the corner is the global node id.
  If the node with the id is not in the triangle, do nothing.

  \param id the global id of the node to be moved
  \param d the amount of movement
*/
void Triangle::moveNode(int id,const Vec3& d)
{
  // move node
  if(id==m_id0) { // move m_p0
    m_p0+=d;
    m_p1-=d;
    m_p2-=d;
  } else if (id==m_id1){ // move p1
    m_p1+=d;
  } else if (id==m_id2){ // move p2
    m_p2+=d;
  } else {
    std::cerr << "trying to move node not in triangle!" << std::endl;
  }
  // recalculate normal and invtrans
  try{
    m_normal=cross(m_p2,m_p1).unit_s();
    m_trans=Matrix3(m_p1,m_p2,m_normal);
    m_invtrans=m_trans.inv();
  }
  catch(VecErr V){ // if unit_s thows exception -> v1||v2
    throw TriangleError();
  }
}

/*!
  Move (translate) whole triangle. 

  \param d the amount of movement
*/
void Triangle::move(const Vec3& d)
{
  m_p0+=d;
}

/*!
	rotate the whole triangle around an axis

    \param origin a point on the axis
    \param axis the orientation of the rotation axis
    \param angle the rotation angle
*/
void Triangle::rotate(const Vec3& origin, const Vec3& axis, double angle)
{
    // get positions of corners 1 & 2
    Vec3 c1=m_p0+m_p1;
    Vec3 c2=m_p0+m_p2;
    // rotate points
    m_p0.rotateBy(origin,axis,angle);
    c1.rotateBy(origin,axis,angle);
    c2.rotateBy(origin,axis,angle);
    // reconstruct edge vectors
    m_p1=c1-m_p0;
    m_p2=c2-m_p0;
    // recalculate normal and invtrans
    try{
        m_normal=cross(m_p2,m_p1).unit_s();
        m_trans=Matrix3(m_p1,m_p2,m_normal);
        m_invtrans=m_trans.inv();
    }
    catch(VecErr V){ // if unit_s thows exception -> v1||v2
        throw TriangleError();
    }
}


/*!
  Transform a point in local coordinates into global coordiantes.
  The local coordinate systems is formed by (P1-P0,P2-P0,N).

  \param p the point to be transformed 
*/
Vec3 Triangle::toGlobal(const Vec3& p)
{
  return m_p0+m_trans*p;
}

/*!
  Transform a point in global coordinates into local coordiantes.
  The local coordinate systems is formed by (P1-P0,P2-P0,N).

  \param p the point to be transformed 
*/
Vec3 Triangle::toLocal(const Vec3& p)
{
  return (p-m_p0)*m_invtrans;
}

/*!
  Get pressure on the triangle from interaction forces
*/
double Triangle::getPressure() const
{
  // calculate normal force
  double F_n=m_force*m_normal;
  // calculate area
  double A=0.5*cross(m_p1,m_p2).norm();
  // retrun force per area
  return F_n/A;  
}
 
/*!
  Get the triangle member function which returns a vector field of a given name. Returns
  NULL if a field with that name doesn't exist.

  \param name the name of the field 
*/
Triangle::VectorFieldFunction Triangle::getVectorFieldFunction(const string& name)
{
  Triangle::VectorFieldFunction f;

  if(name=="force"){
    f=&Triangle::getForce;
  } else {
    f=NULL;
    cerr << "ERROR - invalid name for triangle vector access function" << endl; 
  }

  return f;
}

/*!
  Get the triangle member function which returns a scalar field of a given name. Returns
  NULL if a field with that name doesn't exist.

  \param name the name of the field 
*/
Triangle::ScalarFieldFunction Triangle::getScalarFieldFunction(const string& name)
{
  Triangle::ScalarFieldFunction f;

  if(name=="pressure"){
    f=&Triangle::getPressure;
  } else {
    f=NULL;
    cerr << "ERROR - invalid name for triangle scalar access function" << endl; 
  }

  return f;
}

/*!
  output Triangle to ostream
*/
ostream& operator<<(ostream& ost,const Triangle& T)
{
  ost << "Triangle: (" << T.m_p0 << ") - (" << T.m_p0+T.m_p1 << ") - (" << T.m_p0+T.m_p2 << ") Normal: (" << T.m_normal << ")"; 
  return ost;
}
