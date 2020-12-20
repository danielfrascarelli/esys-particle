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

#ifndef __ROTFRICTIONINTERACTION_H
#define __ROTFRICTIONINTERACTION_H

// -- project includes --
#include "Model/RotPairInteraction.h"
#include "Model/RotParticle.h"
#include "Model/IGParam.h"
#include "Foundation/vec3.h"

// -- I/O includes --
#include <iostream>
using std::ostream;


//double calc_angle( double , double ) ;


/*!
  \struct CRotFrictionIGP
  \brief Interaction parameters for frictional interaction between rotational particles
  \author Shane Latham, Steffen Abe
  $Revision$
  $Date$
*/
class CRotFrictionIGP : public AIGParam
{
public:
  CRotFrictionIGP();
  
  CRotFrictionIGP(
    const  std::string &name,
    double k,
    double mu_d,
    double mu_s,
    double k_s,
    double dt,
    bool   scaling,
    bool   rigid,
    bool   meanR_scaling
  );

  CRotFrictionIGP(
    const  std::string &name,
    double youngsModulus,
    double poissonsRatio,
    double mu_d,
    double mu_s,
    double dt,
    bool   rigid,
    bool   meanR_scaling
  );

  virtual std::string getTypeString() const
  {
    return "RotFriction";
  }

  void setTimeStepSize(double dt);

  double k;
  double mu_d;    // sliding frictional coefficient
  double mu_s;     // max static frictional coefficient
  double k_s;
  double dt;
  bool   scaling;
  bool   rigid;
  bool   meanR_scaling;
};

/*!
   \class CRotFrictionInteraction
   \brief Frictional+Elastic interaction between particles between rotational particles
   \author Shane Latham, Steffen Abe
   $Revision$
   $Date$
*/
class CRotFrictionInteraction : public ARotPairInteraction
{
 public: // types
  typedef CRotFrictionIGP ParameterType;

  typedef double (CRotFrictionInteraction::* ScalarFieldFunction)() const;
  typedef pair<bool,double> (CRotFrictionInteraction::* CheckedScalarFieldFunction)() const;  
  typedef Vec3 (CRotFrictionInteraction::* VectorFieldFunction)() const;

  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);
  
  // type & dummy implementation for parameter setting function 
  typedef void (CRotFrictionInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

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
  bool m_scaling; //!< toggles scaling of elastic properties by particle size
  bool m_rigid; //!< toggles whether to use rigid body friction interactions
  bool m_meanR_scaling; //!< toggles whether to use the mean particle radius or minimum particle radius to define bond radius 
  
  //Quaternion m_init_q1, m_init_q2;
  //Vec3 m_init_pos1 , m_init_pos2;

 public:
  CRotFrictionInteraction();
  CRotFrictionInteraction(CRotParticle*,CRotParticle*,const CRotFrictionIGP&);
  virtual ~CRotFrictionInteraction();

  static string getType() {return "RotFriction";};

  virtual void calcForces();
  virtual void calcSimpleForces();
  virtual void calcRigidBodyForces();
  virtual bool isPersistent();

  void setTimeStepSize(double dt);
  
  void calcNormalForce();
  double getAbsForceDeficit()const;
  double getPotentialEnergy()const;
  double getSlipping()const;
  double getSticking()const;
  double getDissipatedEnergy() const;
  double getAbsSlip() const;
  virtual double Count() const;
  virtual Vec3 getPos() const {return m_cpos;};
  Vec3 getForce() const;
  Vec3 getNormalForce() const;
  Vec3 getTangentialForce() const;

  friend ostream& operator<<(ostream&,const CRotFrictionInteraction&);
  friend class TML_PackedMessageInterface;

  // checkpointing 
  virtual void saveRestartData(std::ostream &oStream);
  virtual void loadRestartData(std::istream &iStream);
};

#endif //__ROTFRICTIONINTERACTION_H
