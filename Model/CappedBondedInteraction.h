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

#ifndef __CAPPEDBONDEDINTERACTION_H
#define __CAPPEDBONDEDINTERACTION_H

// -- project includes --
#include "Model/IGParam.h" // keep this one first - it drags in mpi.h
#include "Model/BondedInteraction.h"
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
  \brief Interaction parameters for bonded interaction with a force limit
  \author Steffen Abe

  $Revision: 894 $
  $Date: 2006-01-19 10:58:58 +0000 (Thu, 19 Jan 2006) $
*/
class CCappedBondedIGP : public CBondedIGP 
{

public:
  double m_force_limit; // maximum force

  CCappedBondedIGP();
  CCappedBondedIGP(const std::string &name, int tag, double normalK, double breakDistance,double forceLimit);

  virtual std::string getTypeString() const  {return "CappedBonded"; }
};

/*!
   \brief Elastic interaction with force limit between bonded particles 
   \author Steffen Abe

   $Revision: 894 $
   $Date: 2006-01-19 10:58:58 +0000 (Thu, 19 Jan 2006) $
*/
class CCappedBondedInteraction : public CBondedInteraction
{
 public: // types
  typedef CCappedBondedIGP ParameterType;

  typedef double (CCappedBondedInteraction::* ScalarFieldFunction)() const;
  typedef pair<bool,double> (CCappedBondedInteraction::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (CCappedBondedInteraction::* VectorFieldFunction)() const;
  
  // type & dummy implementation for parameter setting function 
  typedef void (CCappedBondedInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

 protected:
  double m_force_limit; //!< maximum  allowed force   
  CCappedBondedInteraction(CParticle*,CParticle*);

 public:
  CCappedBondedInteraction();
  CCappedBondedInteraction(
    CParticle *particle1,
    CParticle *particle2,
    const CCappedBondedIGP &params
  );

  virtual ~CCappedBondedInteraction();

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  static string getType() {return "CappedBonded";};

  virtual void calcForces();

  friend class TML_PackedMessageInterface;
};

#endif //__CAPPEDBONDEDINTERACTION_H
