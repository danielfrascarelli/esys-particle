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

#ifndef __WALL_H
#define __WALL_H

//--- project includes ---
#include "Fields/VectorWallFieldSlave.h"
#include "Foundation/vec3.h"
#include "Foundation/console.h"


//--- IO includes ---
#include <iostream>

using std::ostream;
using std::endl;
using std::flush;

class TML_comm;

/*!
  \class CWall
  \brief base class for all walls

  \author Steffen Abe
  $Revision$
  $Date$  
*/
class CWall 
{
protected:
  Vec3 m_origin,m_normal;
  Vec3 m_force;
  Vec3 m_oldpos;
  Vec3 m_vel;
  
public:
  typedef Vec3 (CWall::* VectorFieldFunction)() const;

  CWall();
  CWall(const Vec3&,const Vec3&);
  virtual ~CWall(){};

  void moveBy(const Vec3& v)
  {
    console.XDebug() << "CWall::moveBy: v = " << v << "\n";
    console.XDebug() << "CWall::moveBy: oldpos  = " << m_oldpos << "\n";  
    console.XDebug() << "CWall::moveBy: pre move origin  = " << m_origin << "\n";  
    m_origin += v;
    console.XDebug() << "CWall::moveBy: post move origin = " << m_origin << "\n";      
  };
  void moveTo(const Vec3& v){m_origin=v;};
  void setNormal(const Vec3& v){m_normal=v;};
  void setVel(const Vec3& v){m_vel=v;};
  Vec3 getVel(){return m_vel;};
  inline const Vec3& getOrigin()const {return m_origin;};
  inline const Vec3& getNormal()const {return m_normal;};
  inline void addForce(const Vec3& force){m_force-=force;};
  inline void zeroForce(){m_force=Vec3(0.0,0.0,0.0);};
  inline const Vec3& getForce(){return m_force;};
  inline const Vec3& getPos(){return m_origin;};
  Vec3 getPos() const {return m_origin;};
  Vec3 getForce() const {return m_force;};
  inline double getDisplacement(){return (m_origin-m_oldpos).norm();};
  inline Vec3 getTotalDisplacement(){return (m_origin-m_oldpos);};
  inline void resetDisplacement(){m_oldpos=m_origin;};

  static VectorFieldFunction getVectorFieldFunction(const string&);
  VectorWallFieldSlave<CWall>* generateVectorFieldSlave(TML_Comm*,const string&);
  int getFieldSummationFlag(const string&);

  virtual void writeCheckPoint(ostream&,const string&) const;
  virtual void loadCheckPoint(istream&);

  friend ostream& operator<<(ostream&,const CWall&);
};

#endif //__WALL_H
