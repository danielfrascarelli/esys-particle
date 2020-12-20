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

#ifndef MODEL_ROTLOCALDAMPING_H
#define MODEL_ROTLOCALDAMPING_H

// -- project includes --
#include "Model/LocalDampingIGP.h"
#include "Foundation/vec3.h"
#include "Foundation/quintuple.h"

class CVarMPIBuffer;
class AMPIBuffer;

/*!
  \brief Local rotational damping of the particle motion by a damping coefficient
*/

template <class T>
class CRotLocalDamping
{
protected:
  T *m_p; //!< the particle
  double m_visc; //!< damping coefficient
  double m_dt; //!< time step
  double m_E_diss; //!< dissipated energy
  Vec3 m_force;   //!< current force

public:
  typedef CLocalDampingIGP ParameterType;

  typedef double (CRotLocalDamping::* ScalarFieldFunction)() const;
  typedef pair<bool,double> (CRotLocalDamping::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (CRotLocalDamping::* VectorFieldFunction)() const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);
  
  // type & dummy implementation for parameter setting function 
  typedef void (CRotLocalDamping<T>::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

  CRotLocalDamping(T*,double,double); // to be obsoleted
  CRotLocalDamping(T*,const CLocalDampingIGP&);
  CRotLocalDamping(T*,CLocalDampingIGP*);
  virtual ~CRotLocalDamping();

  void setTimeStepSize(double dt);
  virtual void calcForces();

  virtual bool hasTag(int,int) const;
  virtual Vec3 getPosFirst() const {return m_p->getPos();};
  virtual Vec3 getPosSecond() const {return Vec3(0.0,0.0,0.0);};
  virtual Vec3 getPos() const {return m_p->getPos();};
  vector<int> getAllID() const;

  esys::lsm::quintuple<Vec3,double,Vec3,double,Vec3> getRaw2Data() const
  {
    return 
      esys::lsm::quintuple<Vec3,double,Vec3,double,Vec3>(
        m_p->getPos(),
        m_p->getRad(),
        Vec3::ZERO,
        0,
        getPos()
      );
  }

  double getDissipatedEnergy() const;
  Vec3   getForce() const;
};

#include "Model/RotLocalDamping.hpp"

#endif //__ROTLOCALDAMPING_H
