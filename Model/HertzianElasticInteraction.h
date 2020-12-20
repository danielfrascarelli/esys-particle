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

#ifndef __HERTZIANELASTICINTERACTION_H
#define __HERTZIANELASTICINTERACTION_H

#include "Model/IGParam.h"
#include "Model/Interaction.h"
#include "Model/Particle.h"


/*!
  \brief Interaction group parameters for Hertzian elastic interactions
*/
class CHertzianElasticIGP : public AIGParam
{
public:
  double m_E;
  double m_nu; // poisson ratio

  CHertzianElasticIGP();
  CHertzianElasticIGP(const std::string&,double,double);
  
  virtual std::string getTypeString() const {return "HertzianElastic";}
};

/*!
  \class CHertzianElasticInteraction
  \brief Hertzian Elastic Interaction between free particles
  \author Steffen Abe
  $Revision: 1201 $
  $Date: 2009-07-31 12:25:45 +0200 (Fri, 31 Jul 2009) $
*/
class CHertzianElasticInteraction : public APairInteraction
{
public:

  typedef double (CHertzianElasticInteraction::* ScalarFieldFunction)() const;
  typedef Vec3 (CHertzianElasticInteraction::* VectorFieldFunction)() const;
  typedef pair<bool,double> (CHertzianElasticInteraction::* CheckedScalarFieldFunction)() const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  
  // type & dummy implementation for parameter setting function 
  typedef void (CHertzianElasticInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

private:
  double m_E;//!< Young's modulus
  double m_nu; //!< Poisson ratio
  Vec3   m_force; // caching force for E_pot
  double m_dn; // caching displacement for E_pot
  Vec3   m_cpos; // center position

public:
  typedef CHertzianElasticIGP ParameterType;

  CHertzianElasticInteraction(CParticle*,CParticle*,const CHertzianElasticIGP&);
  virtual ~CHertzianElasticInteraction(){};

  virtual Vec3 getPos() const {return m_cpos;};
  double getPotentialEnergy() const;

  virtual void calcForces();
  Vec3 getForce() const;
};
#endif //__HERTZIANELASTICINTERACTION_H
