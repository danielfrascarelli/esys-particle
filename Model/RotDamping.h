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

#ifndef __ROT_DAMPING_H
#define __ROT_DAMPING_H

// -- project includes --
#include "Model/DampingIGP.h"
#include "Foundation/vec3.h"
#include "Foundation/quintuple.h"
using esys::lsm::quintuple;

/*!
  \brief Damping of the rotational part of the particle motion by an artificial
  viscosity
*/
template <class T>
class CRotDamping
{
protected:
  T *m_p; //!< the particle
  Vec3 m_vref; //!< reference velocity
  double m_visc; //!< artificial viscosity
  double m_dt;   //!< time step
  int m_maxiter;   //!< iteration limit
  double m_E_diss; //!< dissipaed energy
  Vec3 m_force;   //!< current force

  static double s_limit2; //!< square error limit for iteration 
  static int s_flops;       

public:
  typedef CDampingIGP ParameterType;

  typedef double (CRotDamping::* ScalarFieldFunction)() const;
  typedef pair<bool,double> (CRotDamping::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (CRotDamping::* VectorFieldFunction)() const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);
  
  // type & dummy implementation for parameter setting function 
  typedef void (CRotDamping<T>::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

  CRotDamping(T*,CDampingIGP*);
  virtual ~CRotDamping();

  inline void setLimit(double limit){s_limit2=limit*limit;};
  virtual void calcForces();
  void setTimeStepSize(double dt);
  virtual bool hasTag(int,int) const;
  virtual Vec3 getPosFirst() const {return m_p->getPos();};
  virtual Vec3 getPosSecond() const {return Vec3(0.0,0.0,0.0);};
  virtual Vec3 getPos() const {return m_p->getPos();};
  vector<int> getAllID() const;
  quintuple<Vec3,double,Vec3,double,Vec3> getRaw2Data() const
  {
    return 
      quintuple<Vec3,double,Vec3,double,Vec3>(
        m_p->getPos(),
        m_p->getRad(),
        Vec3::ZERO,
        0,
        getPos()
      );
  }

  static void zeroFlops(){s_flops=0;};
  static int Flops(){return s_flops;};

  double getDissipatedEnergy() const;
  Vec3   getForce() const;
};

#include "Model/RotDamping.hpp"

#endif //__ROT_DAMPING_H
