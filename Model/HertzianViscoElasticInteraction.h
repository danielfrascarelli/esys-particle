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

#ifndef __HERTZIANVISCOELASTICINTERACTION_H
#define __HERTZIANVISCOELASTICINTERACTION_H

#include "Model/IGParam.h"
#include "Model/Interaction.h"
#include "Model/Particle.h"


/*!
  \brief Interaction group parameters for Hertzian viscoelastic interactions
*/
class CHertzianViscoElasticIGP : public AIGParam
{
public:
  double m_A;
  double m_E;
  double m_nu; // poisson ratio

  CHertzianViscoElasticIGP();
  CHertzianViscoElasticIGP(const std::string&,double,double,double);
  
  virtual std::string getTypeString() const {return "HertzianViscoElastic";}
};

/*!
  \class CHertzianViscoElasticInteraction
  \brief Hertzian ViscoElastic Interaction between free particles
  \author Laura Heredia & Pablo Richeri
  $Revision: 1 $
  $Date: 2009-12-13 19:00:00 -0300 (Sun, 13 Dec 2009) $
*/
class CHertzianViscoElasticInteraction : public APairInteraction
{
public:

  typedef
    double (CHertzianViscoElasticInteraction::* ScalarFieldFunction)() const;
  typedef
    Vec3 (CHertzianViscoElasticInteraction::* VectorFieldFunction)() const;
  typedef 
    pair<bool,double>
    (CHertzianViscoElasticInteraction::* CheckedScalarFieldFunction)() const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(
    const string&
  );

  // type & dummy implementation for parameter setting function 
  typedef void (CHertzianViscoElasticInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};


private:
  double m_A;//!< Damping constant
  double m_E;//!< Young's modulus
  double m_nu; //!< Poisson ratio
  Vec3   m_force; // caching force for E_pot
  double m_dn; // caching displacement for E_pot
  Vec3   m_cpos; // center position

public:
  typedef CHertzianViscoElasticIGP ParameterType;

  CHertzianViscoElasticInteraction(
    CParticle*,
    CParticle*,
    const CHertzianViscoElasticIGP&
  );
  virtual ~CHertzianViscoElasticInteraction(){};

  virtual Vec3 getPos() const {return m_cpos;};
  double getPotentialEnergy() const;

  virtual void calcForces();
  Vec3 getForce() const;
};
#endif //__HERTZIANVISCOELASTICINTERACTION_H
