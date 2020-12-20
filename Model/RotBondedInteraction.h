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

#ifndef __ROTBONDEDINTERACTION_H
#define __ROTBONDEDINTERACTION_H

// -- project includes --
#include "Model/ARotBondedInteraction.h"
#include "Model/RotParticle.h"
#include "Model/BondedInteractionCpData.h"
#include "Foundation/vec3.h"

// -- I/O includes --
#include <iostream>
using std::ostream;

/*!
  \struct CRotBondedIGP
  \brief Interaction parameters for bonded interaction between rotational particles
  \author Shane Latham, Steffen Abe
  $Revision$
  $Date$
*/

class CRotBondedIGP : public ARotBondedIGP
{
public:
  CRotBondedIGP();
  CRotBondedIGP(
    const  std::string &name,
    double kr,
    double ks,
    double kt,
    double kb,
    double max_nForce,
    double max_shForce,
    double max_tMoment,
    double max_bMoment,
    int    tag,
    bool   scaling,
    bool   AmeanR_scaling,
    double truncated,
    double beta1,
    double beta2
  );

  CRotBondedIGP(
    const  std::string &name,
    double youngsModulus,
    double poissonsRatio,
    double cohesion,
    double tanAngle,
    int    tag,
    bool   AmeanR_scaling,
    double truncated,
    double beta1,
    double beta2
  );

  virtual std::string getTypeString() const
  {
    return "RotBonded";
  }
  
  double max_nForce, max_shForce,max_tMoment, max_bMoment;
  bool   meanR_scaling;
  double truncated;
  double beta1;
  double beta2;
};

/*!
   \class CRotBondedInteraction
   \brief Elastic interaction between bonded particles between rotational particles
   \author Shane Latham, Steffen Abe
   $Revision$
   $Date$
*/
class CRotBondedInteraction : public ARotBondedInteraction
{
 public: // types
  typedef CRotBondedIGP ParameterType;

  /**
   * Used by PIS to save/load check-point data for objects of this type.
   */
  typedef BondedInteractionCpData CheckPointable;

  typedef double (CRotBondedInteraction::* ScalarFieldFunction)() const; 
  typedef pair<bool,double> (CRotBondedInteraction::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (CRotBondedInteraction::* VectorFieldFunction)() const; 

  // type & dummy implementation for parameter setting function 
  typedef void (CRotBondedInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

 protected:
    // maximum force / moment paramaters for failure criteria
  double m_max_nForce;   // always >0
  double m_max_shForce ;
  double m_max_tMoment ;
  double m_max_bMoment ;

  bool m_meanR_scaling; // scale vs. Rmean instead of Rmin
  
 // parameters for truncated M-C criterion
  double m_truncated; 
  double m_beta1;
  double m_beta2;

 public:

  CRotBondedInteraction();
  CRotBondedInteraction(CRotParticle*,CRotParticle*,const CRotBondedIGP&);
  virtual ~CRotBondedInteraction();
                                                                                
  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  static string getType(){return "RotBonded";};

  virtual bool broken();
  double getCriterion() const;
  
  friend ostream& operator<<(ostream&,const CRotBondedInteraction&);
  friend class TML_PackedMessageInterface;
  
  virtual void saveCheckPointData(std::ostream &oStream);
  virtual void loadCheckPointData(std::istream &iStream);

  // save/load of restart parameters
  virtual void saveRestartData(std::ostream &oStream);
  virtual void loadRestartData(std::istream &iStream);
};

#endif //__BONDEDINTERACTION_H
