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

#ifndef __SHORTBONDEDINTERACTION_H
#define __SHORTBONDEDINTERACTION_H

// -- project includes --
#include "Model/IGParam.h" // keep this one first - it drags in mpi.h
#include "Model/Interaction.h"
#include "Model/Particle.h"
#include "Model/BondedInteraction.h" 
#include "Model/ShortBondedInteractionCpData.h"
#include "Foundation/vec3.h"


/*!
  \brief class for a "short" bonded interaction
  
  A bonded interaction where the equilibrium distance is not determined by
  the radii of the particles but by the initial distance, 
  i.e. it allows for overlapping particles. 
  Uses the same parameter class as "normal" bonded interactions
*/
class CShortBondedInteraction : public CBondedInteraction
{
 public: // types
  typedef CBondedIGP ParameterType;
  /**
   * Used by PIS to save/load check-point data for objects of this type.
   */
  typedef ShortBondedInteractionCpData CheckPointable;

  typedef double (CShortBondedInteraction::* ScalarFieldFunction)() const;
  typedef pair<bool,double> (CShortBondedInteraction::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (CShortBondedInteraction::* VectorFieldFunction)() const;
  
  // type & dummy implementation for parameter setting function 
  typedef void (CShortBondedInteraction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

 private:
  
 public:
  CShortBondedInteraction();
  CShortBondedInteraction(CParticle*,CParticle*,const CBondedIGP&);

  virtual ~CShortBondedInteraction();

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  static string getType() {return "ShortBonded";};
  double getEquiDist() const {return m_r0;};

  void saveCheckPointData(std::ostream &oStream);
  void loadCheckPointData(std::istream &iStream);

  friend ostream& operator<<(ostream&,const CBondedInteraction&);
  friend class TML_PackedMessageInterface;
};


#endif // __SHORTBONDEDINTERACTION_H
