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

#ifndef __ROTTHERMBONDEDINTERACTION_H
#define __ROTTHERMBONDEDINTERACTION_H

// -- project includes --
#include "Model/RotThermPairInteraction.h"
#include "Model/RotThermParticle.h"
#include "Model/InteractionParam.h"
#include "Model/BondedInteractionCpData.h"
#include "Foundation/vec3.h"
#include "Model/IGParam.h"

// -- I/O includes --
#include <iostream>
using std::ostream;

/*!
  Interaction parameters for bonded interaction between thermal, rotational particles
*/

 double calc_angle( double , double ) ;

struct CRotThermBondedIGP : public AIGParam
{
  CRotThermBondedIGP();
  CRotThermBondedIGP(
    const std::string &name,
    double kr,
    double ks,
    double kt,
    double kb,
    double max_nForce,
    double max_shForce,
    double max_tMoment,
    double max_bMoment,
    double diffusivity,
    int tag
  );

  double kr,ks,kt,kb ;
  double max_nForce, max_shForce,max_tMoment, max_bMoment;
  double diffusivity ;
  int tag;

  virtual std::string getTypeString() const
  {
    return "RotThermBonded";
  }
};

/*!
   Interaction between bonded, thermal,
   rotational particles
*/
class CRotThermBondedInteraction : public ARotThermPairInteraction
{
 public: // types
  typedef CRotThermBondedIGP ParameterType;

  /**
   * Used by PIS to save/load check-point data for objects of this type.
   */
  typedef BondedInteractionCpData CheckPointable;

  typedef double (CRotThermBondedInteraction::* ScalarFieldFunction)() const; 
  typedef pair<bool,double> (CRotThermBondedInteraction::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (CRotThermBondedInteraction::* VectorFieldFunction)() const; 
  
  // type & dummy implementation for parameter setting function 
  typedef void (CRotThermBondedInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

 private:

  //   protected:
  double m_dist;  //!< current distance, cached from last calcForces()
  double m_min_r; 
  double m_kr ;     //!< spring constant
  double m_ks ;
  double m_kb ;
  double m_kt ;
  
  double m_diffusivity ;  


  double m_max_nForce;   // always >0
  double m_max_shForce ;
  double m_max_tMoment ;
  double m_max_bMoment ;

  double m_nForce;  //  >0, pulling; <0 , compressing
  double m_shForce ; // always >0
  double m_tMoment ;
  double m_bMoment ;

  Vec3 m_force;   //!< current force, cached for E_pot calculation
  Vec3 m_moment ;

  Vec3 m_cpos; // ?
  int m_tag;

 public:

  CRotThermBondedInteraction();
  CRotThermBondedInteraction(CRotThermParticle*,CRotThermParticle*,const CRotThermBondedIGP&);
  virtual ~CRotThermBondedInteraction();
                                                                                
  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  static string getType(){return "RotThermBonded";};

  int getTag() const;
  void setTag(int tag);

  void calcForces();
  void calcHeatTrans() ;
  //void setBreak(double);
  bool broken();

//26/Mar added !
  Vec3 getBondedVector1() const;
  Vec3 getBondedVector2() const;

  double getPotentialEnergy() const;
  double getNormalPotentialEnergy() const;
  double getShearPotentialEnergy() const;
  double getTwistPotentialEnergy() const;
  double getBendPotentialEnergy() const;
  double getCriterion() const;
  Vec3   getForce() const;
  virtual Vec3 getPos() const {return m_cpos;};

  Vec3 getCentrePtDiff() const;
  Vec3 getInitialCentrePtDiff() const;
  Vec3 getInitialMidPoint() const;

  Vec3 getShearDiff() const;

  friend ostream& operator<<(ostream&,const CRotThermBondedInteraction&);
  friend class TML_PackedMessageInterface;
  
  virtual void saveCheckPointData(std::ostream &oStream);

  virtual void loadCheckPointData(std::istream &iStream);

  // save/load of restart parameters
  virtual void saveRestartData(std::ostream &oStream);
  virtual void loadRestartData(std::istream &iStream);
};

#endif //__BONDEDINTERACTION_H
