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

#ifndef __ROTTHERMFRICTIONINTERACTION_H
#define __ROTTHERMFRICTIONINTERACTION_H

// -- project includes --
#include "Model/RotThermPairInteraction.h"
#include "Model/RotThermParticle.h"
#include "Model/InteractionParam.h"
#include "Foundation/vec3.h"
#include "Model/IGParam.h"

// -- I/O includes --
#include <iostream>
using std::ostream;


//double calc_angle( double , double ) ;


/*!
  Interaction parameters for frictional interaction between Thermal,rotational particles
*/
struct CRotThermFrictionIGP : public AIGParam
{
  CRotThermFrictionIGP();

  CRotThermFrictionIGP(
    const std::string &name,
    double k,
    double mu_d,
    double mu_s,
    double k_s,
    double diffusivity,
    double dt
  );

  double k;
  double mu_d;    // sliding frictional coefficient
  double mu_s;     // max static frictional coefficient
  double k_s;
  double dt;
  double diffusivity ;

  virtual std::string getTypeString() const
  {
    return "RotThermFriction";
  }

  inline void setTimeStepSize(double deltaT)
  {
    dt = deltaT;
  }
};

/*!
   Frictional+Elastic interaction between particles between thermal ,rotational particles
*/
class CRotThermFrictionInteraction : public ARotThermPairInteraction
{
 public: // types
  typedef CRotThermFrictionIGP ParameterType;

  typedef double (CRotThermFrictionInteraction::* ScalarFieldFunction)() const;
  typedef pair<bool,double> (CRotThermFrictionInteraction::* CheckedScalarFieldFunction)() const;  
  typedef Vec3 (CRotThermFrictionInteraction::* VectorFieldFunction)() const;

  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);
  
  // type & dummy implementation for parameter setting function 
  typedef void (CRotThermFrictionInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

  inline void setTimeStepSize(double deltaT)
  {
    m_dt = deltaT;
  }
// protected:
 private:

  double m_k;     //!< spring constant
  double m_r0;    //!< equilibrium distance
  double m_mu_d;  //!< coefficient of dynamic friction
  double m_mu_s;  //!< coefficient of static friction
  double m_ks;    //!< shear stiffness (Cundall)
  double m_dt;    //!< time step
  Vec3 m_Ffric;   //!< current frictional force
  Vec3 m_force_deficit; //!< difference between fric. force & force necessary for slip
  Vec3 m_cpos; //!< contact position
  Vec3 m_normal_force; //!< current normal force
  bool m_is_slipping; //!< static/dynamic status of the interaction
  bool m_is_touching; //!< contact status of the interaction
  double m_E_diss; //!< dissipated energy
  double m_diffusivity; //!< thermal diffusivity
  Vec3 m_ds; //!< tangitial displacement at this time step
  //Quaternion m_init_q1, m_init_q2;
  //Vec3 m_init_pos1 , m_init_pos2;

 public:
  CRotThermFrictionInteraction();
  CRotThermFrictionInteraction(CRotThermParticle*,CRotThermParticle*,const CRotThermFrictionIGP&);
  virtual ~CRotThermFrictionInteraction();

  static string getType() {return "RotThermFriction";};

  virtual void calcForces();
  void calcHeatFrict();
  void calcHeatTrans();
  virtual bool isPersistent();

  void calcNormalForce();
  double getAbsForceDeficit()const;
  double getPotentialEnergy()const;
  double getSlipping()const;
  double getTouching()const;
  double getSticking()const;
  double getDissipatedEnergy() const;
  inline Vec3 getDs() {return m_ds;}
  virtual Vec3 getPos() const {return m_cpos;}
  Vec3 getForce() const;
  Vec3 getNormalForce() const;

  friend ostream& operator<<(ostream&,const CRotThermFrictionInteraction&);
  friend class TML_PackedMessageInterface;

  // checkpointing 
  virtual void saveRestartData(std::ostream &oStream);
  virtual void loadRestartData(std::istream &iStream);
};

#endif //__ROTFRICTIONINTERACTION_H
