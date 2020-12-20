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

#ifndef __ROTELASTICINTERACTION_H
#define __ROTELASTICINTERACTION_H

#include "Model/RotPairInteraction.h"
#include "Model/RotParticle.h"
#include "Model/IGParam.h"


/*!
  \brief Interaction group parameters for CRotElasticInteractionGroups
*/
class CRotElasticIGP : public AIGParam
{
 protected:
 public:
  CRotElasticIGP();
  CRotElasticIGP(const std::string &name, double kr, bool scaling);
  
  double m_kr;
  bool   m_scaling;

  virtual void packInto(CVarMPIBuffer*) const;
  void setNormalSpringConst(double k){m_kr=k;};
  double getNormalSpringConst() const{return m_kr;};
  
  virtual std::string getTypeString() const {return "RotElastic";}

  friend ostream& operator<<(ostream&,const CRotElasticIGP&);
};

CRotElasticIGP* extractRotElasticIGP(AMPIBuffer*);
CRotElasticIGP* extractRotElasticIGP_p(AMPIBuffer*);

/*!
  \class CRotElasticInteraction
  \brief Elastic Interaction between free rotational particles
*/
class CRotElasticInteraction : public ARotPairInteraction
{
public:
  typedef double (CRotElasticInteraction::* ScalarFieldFunction)() const;
  typedef pair<bool,double> (CRotElasticInteraction::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (CRotElasticInteraction::* VectorFieldFunction)() const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);
  
  // type & dummy implementation for parameter setting function 
  typedef void (CRotElasticInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

private:
  double m_kr; //!< spring constant
  Vec3   m_force; //!< caching force for E_pot
  double m_nForce; //!< normal force, always >= 0
  Vec3   m_cpos; //!< current position
  Vec3   m_D; //!< initial positions of the particles
  bool   m_scaling; //!< scaling of normal stiffness with particle size 

public:
  typedef CRotElasticIGP ParameterType;

  CRotElasticInteraction();
  CRotElasticInteraction(CRotParticle*,CRotParticle*,const CRotElasticIGP&);
  virtual ~CRotElasticInteraction(){};

  virtual Vec3 getPos() const {return m_cpos;}

  static string getType(){return "RotElastic";}

  virtual void calcForces();

  Vec3   getForce() const;
  double getPotentialEnergy() const;

  friend ostream& operator<<(ostream&,const CRotElasticInteraction&);
  
  // save/load of restart parameters
  virtual void saveRestartData(std::ostream &oStream);
  virtual void loadRestartData(std::istream &iStream);
};
#endif //__ELASTICINTERACTION_H
