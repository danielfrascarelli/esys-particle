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

#ifndef MODEL_ABCDAMPING_H
#define MODEL_ABCDAMPING_H

// -- project includes --
#include "Model/Damping.h"
#include "Model/ABCDampingIGP.h"

/*!
  Damping for absorbing boundary conditions - damping increases exponentially
  towards a given plane (boundary)
*/
 template <class ParticleType>
 class ABCDamping : public CDamping<ParticleType>
{
 protected:
  Vec3 m_pos;
  Vec3 m_normal;
  double m_c1;

 public:
  typedef ABCDampingIGP ParameterType;


  typedef double (ABCDamping::* ScalarFieldFunction)() const;
  typedef pair<bool,double> (ABCDamping::* CheckedScalarFieldFunction)() const;
  typedef Vec3 (ABCDamping::* VectorFieldFunction)() const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);
  
  // type & dummy implementation for parameter setting function 
  typedef void (ABCDamping::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

  ABCDamping(ParticleType*,ABCDampingIGP*);
  virtual ~ABCDamping();
};

#include "Model/ABCDamping.hpp"

#endif // MODEL_ABCDAMPING_H
