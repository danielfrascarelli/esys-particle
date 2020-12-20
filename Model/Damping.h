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

#ifndef MODEL_DAMPING_H
#define MODEL_DAMPING_H

// -- project includes --
#include "Model/DampingIGP.h"
#include "Foundation/vec3.h"
#include "Foundation/quintuple.h"

class CVarMPIBuffer;
class AMPIBuffer;

/*!
  \brief Damping of the particle motion by an artificial viscosity
*/

template <class T>
class CDamping
{
protected:
  T *m_p; //!< the particle
  Vec3 m_vref; //!< reference velocity
  double m_visc; //!< artificial viscosity
  double m_dt;   //!< time step
  int m_maxiter;   //!< iteration limit
  double m_E_diss; //!< dissipated energy
  Vec3 m_force;   //!< current force

  static double s_limit2; //!< square error limit for iteration 
  static int s_flops;       

public:
  typedef CDampingIGP ParameterType;

  typedef double (CDamping::* ScalarFieldFunction)() const;
  typedef pair<bool,double> (CDamping::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (CDamping::* VectorFieldFunction)() const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  // type & dummy implementation for parameter setting function 
  typedef void (CDamping::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};
  
  CDamping(T*,const Vec3&,double,double,int); // to be obsoleted
  CDamping(T*,const CDampingIGP&);
  CDamping(T*,CDampingIGP*);
  virtual ~CDamping();

  inline void setLimit(double limit){s_limit2=limit*limit;};
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

  static void zeroFlops(){s_flops=0;};
  static int Flops(){return s_flops;};

  double getDissipatedEnergy() const;
  Vec3   getForce() const;
};

#include "Model/Damping.hpp"

#endif //__DAMPING_H
