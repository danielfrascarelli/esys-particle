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

#ifndef __SPRINGDASHPOTFRICTIONINTERACTION_H
#define __SPRINGDASHPOTFRICTIONINTERACTION_H

#include "Model/IGParam.h"
#include "Model/Interaction.h"
#include "Model/Particle.h"
#include "Foundation/vec3.h"

#include <iostream>
#include <utility>

/*!
  \brief Interaction parameters for frictional interaction
*/
class CSpringDashpotFrictionIGP : public AIGParam
{
public:
  CSpringDashpotFrictionIGP();

  CSpringDashpotFrictionIGP(const std::string &name, double youngModulus, double poissionRatio, double cor, double fricCoef, double dT);

  virtual std::string getTypeString() const {return "Spring_Dashpot_Friction";}

  void setTimeStepSize(double dt);

  double youngModulus;
  double poissonRatio;
  double cor;
  double mu;
  double dt;
};

/*!
   \class CSpringDashpotFrictionInteraction
   \brief Frictional+Elastic+Dashpot interaction between particles
   \author Steffen Abe
   $Revision$
   $Date$
*/
class CSpringDashpotFrictionInteraction : public APairInteraction
{
 public: // types
  typedef CSpringDashpotFrictionIGP ParameterType;

  typedef double (CSpringDashpotFrictionInteraction::* ScalarFieldFunction)() const;
  typedef std::pair<bool,double> (CSpringDashpotFrictionInteraction::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (CSpringDashpotFrictionInteraction::* VectorFieldFunction)() const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);
  
  // type & dummy implementation for parameter setting function 
  typedef void (CSpringDashpotFrictionInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

 protected:
  double m_k;     //!< spring constant
  double m_r0;    //!< equilibrium distance
  double m_mu;    //!< coefficient of friction
  double m_ks;    //!< shear stiffness (Cundall)
  double m_dt;    //!< time step
  double m_damping; //!< damping coefficient

  Vec3 m_Ffric;   //!< current frictional force
  Vec3 m_force_deficit; //!< difference between fric. force & force necessary for slip
  Vec3 m_cpos; //!< contact position
  Vec3 m_normal_force; //!< current normal force
  bool m_is_slipping; //!< static/dynamic status of the interaction
  bool m_is_touching; //!< contact status of the interaction
  double m_E_diss; //!< dissipated energy

 public:
  CSpringDashpotFrictionInteraction();
  CSpringDashpotFrictionInteraction(CParticle*,CParticle*);
  CSpringDashpotFrictionInteraction(CParticle*,CParticle*,const CSpringDashpotFrictionIGP&);
  virtual ~CSpringDashpotFrictionInteraction();

  static string getType() {return "Spring_Dashpot_Friction";};
  
  virtual void calcForces();
  virtual bool isPersistent();

  void setTimeStepSize(double dt);

  std::pair<bool,double> getAbsFrictionalForce() const;
  std::pair<bool,double> getAbsFN() const;
  std::pair<bool,double> getAbsMuFN() const;
  std::pair<bool,double> getSlipVelocity() const;
  std::pair<bool,double> getNormalStress() const;
  std::pair<bool,double> getMaxFricStress() const;
  std::pair<bool,double>  getAbsFrictionalStress() const;

  double getAbsForceDeficit() const;
  double getPotentialEnergy() const;
  double getSlipping()const;
  double getSticking()const;
  double getDissipatedEnergy() const;
  virtual double Count() const;
  Vec3 getForce() const;
  Vec3 getNormalForce() const;
  virtual Vec3 getPos() const {return m_cpos;};

  std::pair<bool,double> getMuEff(const Vec3&,const Vec3&) const;
  std::pair<bool,double> getMuEffXY() const {return getMuEff(Vec3(1.0,0.0,0.0),Vec3(0.0,1.0,0.0));};
  std::pair<bool,double> getMuEffXZ() const {return getMuEff(Vec3(1.0,0.0,0.0),Vec3(0.0,0.0,1.0));};

  friend std::ostream& operator<<(std::ostream&,const CSpringDashpotFrictionInteraction&);
  friend class TML_PackedMessageInterface;
  
  // save/load of restart parameters
  virtual void saveRestartData(std::ostream &oStream);
  virtual void loadRestartData(std::istream &iStream);
};

#endif //__SPRINGDASHPOTFRICTIONINTERACTION_H
