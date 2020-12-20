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

#ifndef __HERTZMINDLINVISCOINTERACTION_H
#define __HERTZMINDLINVISCOINTERACTION_H

#include "Model/IGParam.h"
#include "Model/Interaction.h"
#include "Model/Particle.h"
#include "Foundation/vec3.h"

#include <iostream>
#include <utility>

/*!
  \brief Interaction group parameters for Hertz-Mindlin interactions
*/
class CHertzMindlinViscoIGP : public AIGParam
{
public:
  CHertzMindlinViscoIGP();

  CHertzMindlinViscoIGP(
    const std::string &name,
    double E,
    double nu,
    double fricCoef,
    double restiCoef,
    double dT
  );

  virtual std::string getTypeString() const
  {
    return "HertzMindlinVisco";
  }

  void setTimeStepSize(double dt);

  double m_E;  // Young's modulus
  double m_nu; // Poisson ratio
  double mu;   // Friction coefficient
  double m_COR; //Coefficient of restitution
  double dt;
};

/*!
   \class CHertzMindlinViscoInteraction
   \brief Frictional+ViscoElastic interaction between particles
   \author Qi Shao & Dion Weatherley
   $Date: 2016-12-14 10:00:00 (wed, 14 dec 2016) $
*/
class CHertzMindlinViscoInteraction : public APairInteraction
{
 public: // types
  typedef CHertzMindlinViscoIGP ParameterType;

  typedef
    double (CHertzMindlinViscoInteraction::* ScalarFieldFunction)(
    ) const;
  typedef
    std::pair<bool,double>
    (CHertzMindlinViscoInteraction::* CheckedScalarFieldFunction)(
    ) const;
  typedef
    Vec3 (CHertzMindlinViscoInteraction::* VectorFieldFunction)(
    ) const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(
    const string&
  );
  static VectorFieldFunction getVectorFieldFunction(const string&);
  
  // type & dummy implementation for parameter setting function 
  typedef void (CHertzMindlinViscoInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

 protected:
  double m_E;     //!< Young's modulus
  double m_nu;    //!< Poisson ratio
  double m_r0;    //!< equilibrium distance
  double m_mu;    //!< coefficient of friction
  double m_COR;   //!< coefficient of restitution
  double m_dt;    //!< time step
  Vec3 m_Ffric;   //!< current frictional force
  Vec3 m_force_deficit; //!< difference between fric. force & force necessary for slip
  Vec3 m_cpos; //!< contact position
  Vec3 m_normal_force; //!< current normal force
  double m_dn; //!< caching displacement for E_pot
  bool m_is_slipping; //!< static/dynamic status of the interaction
  bool m_is_touching; //!< contact status of the interaction
  double m_E_diss; //!< dissipated energy
  double m_E_visc_normal; //!< viscous energy increment dissipated in normal direction
  double m_E_visc_shear; //!< viscous energy increment dissipated in shear direction

 public:
  CHertzMindlinViscoInteraction();
  CHertzMindlinViscoInteraction(CParticle*,CParticle*);
  CHertzMindlinViscoInteraction(
    CParticle*,
    CParticle*,
    const CHertzMindlinViscoIGP&
  );
  virtual ~CHertzMindlinViscoInteraction();

  static string getType() {return "CHertzMindlinViscoInteraction";};

  virtual void calcForces();
  virtual bool isPersistent();

  void setTimeStepSize(double dt);

  std::pair<bool,double> getAbsFrictionalForce() const;
  std::pair<bool,double> getAbsFN() const;
  std::pair<bool,double> getAbsMuFN() const;
  std::pair<bool,double> getSlipVelocity() const;
  std::pair<bool,double> getNormalStress() const;
  std::pair<bool,double> getMaxFricStress() const;
  std::pair<bool,double> getAbsFrictionalStress() const;

  double getAbsForceDeficit() const;
  double getPotentialEnergy() const;
  double getSlipping()const;
  double getSticking()const;
  double getDissipatedEnergy() const;
  double getNormalViscousEnergy() const;
  double getShearViscousEnergy() const;
  double Count() const;
  Vec3 getForce() const;
  Vec3 getNormalForce() const;
  virtual Vec3 getPos() const {return m_cpos;};

  std::pair<bool,double> getMuEff(const Vec3&,const Vec3&) const;
  std::pair<bool,double> getMuEffXY() const
  {
    return getMuEff(Vec3(1.0,0.0,0.0),Vec3(0.0,1.0,0.0));
  };
  std::pair<bool,double> getMuEffXZ() const
  {
    return getMuEff(Vec3(1.0,0.0,0.0),Vec3(0.0,0.0,1.0));
  };

  friend std::ostream& operator<<(
    std::ostream&,
    const CHertzMindlinViscoInteraction&
  );
  friend class TML_PackedMessageInterface;
};

#endif //__HERTZMINDLINVISCOINTERACTION_H
