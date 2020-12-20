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

#ifndef __SPHEREBODY_H
#define __SPHEREBODY_H

//--- project includes ---
//#include "Fields/VectorWallFieldSlave.h"
#include "Foundation/vec3.h"
#include "Foundation/console.h"


//--- IO includes ---
#include <iostream>

using std::ostream;
using std::endl;
using std::flush;

class TML_comm;

/*!
  \class CSphereBody
  \brief base class for spherical non-inertial bodies (similar to simple walls)

  \author Dion Weatherley
  $Revision$
  $Date$  
*/
class CSphereBody
{
   protected:
      Vec3 m_centre;
      double m_radius;
      Vec3 m_force;
      Vec3 m_oldpos;
      Vec3 m_vel;

   public:
      CSphereBody();
      CSphereBody(const Vec3&,const double&);
      virtual ~CSphereBody(){};
       
      void moveBy (const Vec3& v) {m_centre += v;};
      void moveTo (const Vec3& v) {m_centre = v;};
      void setVel (const Vec3& v) {m_vel = v;};
      Vec3 getVel () {return m_vel;};

      inline const Vec3& getCentre() const {return m_centre;};
      inline const double& getRadius() const {return m_radius;};
      inline void addForce(const Vec3& force) {m_force -= force;};
      inline void zeroForce(){m_force = Vec3(0.0,0.0,0.0);};
      inline const Vec3& getForce(){return m_force;};
      inline const Vec3& getPos(){return m_centre;};
      //Vec3 getPos(){return m_centre;};
      //ec3 getForce(){return m_force;};
      inline double getDisplacement(){return (m_centre - m_oldpos).norm();};
      inline Vec3 getTotalDisplacement(){return (m_centre - m_oldpos);};
      inline void resetDisplacement(){m_oldpos = m_centre;};

      virtual void writeCheckPoint(ostream&,const string&) const;
      virtual void loadCheckPoint(istream&);

      friend ostream& operator<<(ostream&,const CSphereBody&);
};

#endif //__SPHEREBODY_H
