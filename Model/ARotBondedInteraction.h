/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2019 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0              //
//                                                         //
/////////////////////////////////////////////////////////////

#ifndef __AROTBONDEDINTERACTION_H
#define __AROTBONDEDINTERACTION_H

// -- project includes --
#include "Model/RotPairInteraction.h"
#include "Model/RotParticle.h"
#include "Model/BondedInteractionCpData.h"
#include "Model/IGParam.h"
#include "Foundation/vec3.h"

// -- I/O includes --
#include <iostream>
using std::ostream;


double calc_angle( double , double ) ;

/*!
  \struct ARotBondedIGP
  \brief Interaction parameters for bonded interaction between rotational particles
  \author Shane Latham, Steffen Abe
  
  This structure only contains the parameters common to all rotational bonded interactions.
  Additional parameters for specific interactions are contained in derived classes
  
  $Revision$
  $Date$
*/
class ARotBondedIGP : public AIGParam
{
public:
  ARotBondedIGP();
  ARotBondedIGP(
    const  std::string &name,
    double kr,
    double ks,
    double kt,
    double kb,
    int tag,  
    bool scaling
  );

  ARotBondedIGP(
    const  std::string &name,
    double youngsModulus,
    double poissonsRatio,
    int    tag
  );

  virtual std::string getTypeString() const
  {
    return "ARotBonded";
  }
  
  double kr,ks,kt,kb ;
  int    tag;
  bool   scaling;
};

/*!
   \class ARotBondedInteraction
   \brief Base class for a bonded interaction between bonded particles between rotational particles
   \author Shane Latham, Steffen Abe

   Contains only the part common to all rotational bonded interactions, i.e. force & stress calculations.
   Specifics of parameterisation and bond breaking are handled in derived classes.

   $Revision$
   $Date$
*/
class ARotBondedInteraction : public ARotPairInteraction
{
 public: // types
  typedef ARotBondedIGP ParameterType;

  typedef double (ARotBondedInteraction::* ScalarFieldFunction)() const; 
  typedef pair<bool,double> (ARotBondedInteraction::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (ARotBondedInteraction::* VectorFieldFunction)() const; 

  // type & dummy implementation for parameter setting function 
  typedef void (ARotBondedInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

 protected:

  //   protected:
  double m_dist;  //!< current distance, cached from last calcForces()
  double m_r0;  //!< equilibrium separation

  double m_kr ;     //!< spring constant
  double m_ks ;
  double m_kb ;
  double m_kt ;

  double m_nForce;  //  >0, pulling; <0 , compressing
  double m_shForce ; // always >0
  double m_tMoment ;
  double m_bMoment ;

  Vec3 m_force;   //!< current force, cached for E_pot calculation
  Vec3 m_moment ;

  Vec3 m_cpos; // ?
  Vec3 m_D; //!< initial positions of the particles
  int  m_tag;
  bool m_scaling;
  
 public:

  ARotBondedInteraction();
  ARotBondedInteraction(CRotParticle*,CRotParticle*,const ARotBondedIGP&);
  virtual ~ARotBondedInteraction();
                                                                                
  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  static string getType(){return "RotBonded";};

  int getTag() const;
  void setTag(int tag);

  void calcForces();
  //void setBreak(double);
  virtual bool broken() =0; // pure virtual - defined in derived classes

  double getPotentialEnergy() const;
  double getNormalPotentialEnergy() const;
  double getShearPotentialEnergy() const;
  double getTwistPotentialEnergy() const;
  double getBendPotentialEnergy() const;
  Vec3   getForce() const;
  Vec3   getNormalForce() const;
  Vec3   getTangentialForce() const;
  virtual Vec3 getPos() const {return m_cpos;};

  Vec3 getCentrePtDiff() const;
  Vec3 getInitialCentrePtDiff() const;
  Vec3 getInitialMidPoint() const;

  Vec3 getP2ShearForcePt() const;
  Vec3 getP1ShearForcePt() const;
  Vec3 getContactPoint() const;

  Vec3 getShearDiff() const;

  virtual void saveCheckPointData(std::ostream &oStream)=0;
  virtual void loadCheckPointData(std::istream &iStream)=0;

  // save/load of restart parameters
  virtual void saveRestartData(std::ostream &oStream)=0;
  virtual void loadRestartData(std::istream &iStream)=0;
};

#endif //__BONDEDINTERACTION_H
