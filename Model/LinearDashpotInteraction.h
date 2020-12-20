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

#ifndef __LINEARDASHPOTINTERACTION_H
#define __LINEARDASHPOTINTERACTION_H

#include "Model/IGParam.h"
#include "Model/Interaction.h"
#include "Model/Particle.h"


/*!
  \brief Interaction group parameters for Linear Dashpot interactions
*/
class CLinearDashpotIGP : public AIGParam
{
public:
  double m_damp;
  double m_cutoff;

  CLinearDashpotIGP();
  CLinearDashpotIGP(const std::string&,double,double);
  
  virtual std::string getTypeString() const {return "LinearDashpot";}
};

/*!
  \class CLinearDashpotInteraction
  \brief Linear Dashpot Interaction between free or bonded particles (to be used in addition to an elastic or bonded Interaction, not exclusively)
  \author Steffen Abe
  $Revision: 1201 $
  $Date: 2009-07-31 12:25:45 +0200 (Fri, 31 Jul 2009) $
*/
class CLinearDashpotInteraction : public APairInteraction
{
public:

  typedef double (CLinearDashpotInteraction::* ScalarFieldFunction)() const;
  typedef Vec3 (CLinearDashpotInteraction::* VectorFieldFunction)() const;
  typedef pair<bool,double> (CLinearDashpotInteraction::* CheckedScalarFieldFunction)() const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  
  // type & dummy implementation for parameter setting function 
  typedef void (CLinearDashpotInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

private:
  double m_damp;//!< spring constant
  double m_cutoff; //!< interaction distance cutoff, relative to particle radii
  double m_cross_section; //!< cross section of dashpot, calculated from particle radii
  Vec3   m_force; // caching force for E_pot
  Vec3   m_cpos; // center position

public:
  typedef CLinearDashpotIGP ParameterType;

  CLinearDashpotInteraction(CParticle*,CParticle*,const CLinearDashpotIGP&);
  virtual ~CLinearDashpotInteraction(){};

  virtual Vec3 getPos() const {return m_cpos;};
  double getPotentialEnergy() const;

  virtual void calcForces();
  Vec3 getForce() const;
};
#endif //__HERTZIANELASTICINTERACTION_H
