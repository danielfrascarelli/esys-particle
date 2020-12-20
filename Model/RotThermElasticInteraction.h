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

#ifndef __ROTTHERMELASTICINTERACTION_H
#define __ROTTHERMELASTICINTERACTION_H

#include "Model/RotThermPairInteraction.h"
#include "Model/RotThermParticle.h"
#include "Model/IGParam.h"

/*!
  Interaction group parameters for CRotThermElasticInteractionGroups
*/
class CRotThermElasticIGP : public AIGParam
{
 protected:
 public:
  CRotThermElasticIGP();

  CRotThermElasticIGP(
    const std::string &name,
    double normalK,
    double diffusivity
  );
  
  double m_kr;
  double diffusivity ;

  virtual void  packInto(CVarMPIBuffer*) const;
  void setSpringConst(double k){m_kr=k;};
  double getSpringConst() const{return m_kr;};


  void setDiffusivity(double d){diffusivity = d;};
  double getDiffusivity() const{return diffusivity ; };


  friend ostream& operator<<(ostream&,const CRotThermElasticIGP&);

  virtual std::string getTypeString() const
  {
    return "RotThermElastic";
  }
};

CRotThermElasticIGP* extractRotThermElasticIGP(AMPIBuffer*);
CRotThermElasticIGP* extractRotThermElasticIGP_p(AMPIBuffer*);

/*!
  Elastic Interaction between free thermal, rotational particles
*/
class CRotThermElasticInteraction : public  ARotThermPairInteraction
{
public:
  typedef double (CRotThermElasticInteraction::* ScalarFieldFunction)() const;
  typedef pair<bool,double> (CRotThermElasticInteraction::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (CRotThermElasticInteraction::* VectorFieldFunction)() const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  // type & dummy implementation for parameter setting function 
  typedef void (CRotThermElasticInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

private:
  double m_kr; //!< spring constant
  Vec3   m_force; //!< caching force for E_pot
  double m_nForce; //!< normal force, always >= 0
  Vec3   m_cpos;  //!< current position
  Vec3   m_D; //!< initial positions of the particles
  double m_diffusivity; //!< thermal diffusivity

public:
  typedef CRotThermElasticIGP ParameterType;

  CRotThermElasticInteraction();
  CRotThermElasticInteraction(CRotThermParticle*,CRotThermParticle*,const CRotThermElasticIGP&);
  virtual ~CRotThermElasticInteraction(){};

  virtual Vec3 getPos() const {return m_cpos;};

  static string getType(){return "RotThermElastic";}

  virtual void calcForces();
  void calcHeatTrans() ;
//  void calcHeatFrict(){} ;
     
  Vec3 getForce() const;
  double getPotentialEnergy() const;
//26/Mar added
  Vec3 getBondedVector() const ;

  friend ostream& operator<<(ostream&,const CRotThermElasticInteraction&);
  
  // save/load of restart parameters
  virtual void saveRestartData(std::ostream &oStream);
  virtual void loadRestartData(std::istream &iStream);
};
#endif //__ELASTICINTERACTION_H
