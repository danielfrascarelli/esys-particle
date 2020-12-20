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

#ifndef __BONDEDINTERACTION_H
#define __BONDEDINTERACTION_H

// -- project includes --
#include "Model/IGParam.h" // keep this one first - it drags in mpi.h
#include "Model/Interaction.h"
#include "Model/Particle.h"
#include "Model/BondedInteractionCpData.h"
#include "Foundation/vec3.h"

// -- I/O includes --
#include <iostream>
using std::ostream;

// -- STL includes --
#include <utility>

using std::pair;

/*!
  \brief Interaction parameters for bonded interaction
  \author Steffen Abe

  $Revision$
  $Date$
*/
class CBondedIGP : public AIGParam
{
public:
  CBondedIGP();

  CBondedIGP(const std::string &name, int tag, double normalK, double breakDistance, bool scaling=true);

  virtual std::string getTypeString() const
  {
    return "Bonded";
  }

  double k; //!< Spring constant. 
  double rbreak; //!< Breaking strain
  int tag;
  bool m_scaling;
};

/*!
   \brief Elastic interaction between bonded particles
   \author Steffen Abe

   $Revision$
   $Date$
*/
class CBondedInteraction : public APairInteraction
{
 public: // types
  typedef CBondedIGP ParameterType;
  /**
   * Used by PIS to save/load check-point data for objects of this type.
   */
  typedef BondedInteractionCpData CheckPointable;

  typedef double (CBondedInteraction::* ScalarFieldFunction)() const;
  typedef pair<bool,double> (CBondedInteraction::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (CBondedInteraction::* VectorFieldFunction)() const;
  
  // type & dummy implementation for parameter setting function 
  typedef void (CBondedInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

 protected:
  double m_k;     //!< spring constant
  double m_r0;    //!< equilibrium distance
  double m_dist;  //!< current distance, cached from last calcForces()
  double m_break; //!< breaking distance
  Vec3 m_force;   //!< current force, cached for E_pot calculation
  Vec3 m_cpos;
  int  m_tag;     //!< Interaction tag;
  bool m_scaling; //!< scaling k with particle radius 
  
  CBondedInteraction(CParticle*,CParticle*);

 public:
  CBondedInteraction();
  CBondedInteraction(
    CParticle *particle1,
    CParticle *particle2,
    const CBondedIGP &params
  );

  virtual ~CBondedInteraction();

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  static string getType() {return "Bonded";};

  virtual void calcForces();
  void setBreak(double);
  bool broken();

  inline int getTag() const {return m_tag;}
  inline void setTag(int tag) {m_tag = tag;}

  double getCriterion() const;
  double getPotentialEnergy() const;
  double getStrain() const;
  Vec3 getForce() const;

  virtual Vec3 getPos() const {return m_cpos;};
  virtual void saveCheckPointData(std::ostream &oStream);

  friend ostream& operator<<(ostream&,const CBondedInteraction&);
  friend class TML_PackedMessageInterface;

  // save/load of restart parameters
  virtual void saveRestartData(std::ostream &oStream);
  virtual void loadRestartData(std::istream &iStream);
};

#endif //__BONDEDINTERACTION_H
