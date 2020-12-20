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

#ifndef __FRACTALFRICTION_H
#define __FRACTALFRICTION_H

// -- project includes --
#include "Model/FrictionInteraction.h"

#include <boost/shared_ptr.hpp>

/*!
  \brief Interaction parameters for frictional interaction with a fractal
  distribution of the coefficient of friction
*/
class FractalFrictionIGP : public AIGParam
{
public:
  virtual std::string getTypeString() const {return "FractalFriction";}
  
  void setTimeStepSize(double timeStepSize)
  {
    this->dt = timeStepSize;
  }
  
  double k;
  double mu_0;
  double k_s;
  double dt;
  boost::shared_ptr<double> mu;          //!< pointer to the array of friction coeff.
  double x0,y0,dx,dy;  //!< origin and grid spacing of the array
  int nx,ny;           //!< array size

  FractalFrictionIGP();
  FractalFrictionIGP(const FractalFrictionIGP &);
  ~FractalFrictionIGP();
  
  FractalFrictionIGP &operator=(const FractalFrictionIGP &);
};

/*!
  \brief Frictional+Elastic interaction between particles with fractal distribution of the
  coefficient of friction
*/
class CFractalFriction : public CFrictionInteraction
{
 public: // types
  typedef FractalFrictionIGP ParameterType;

  typedef double (CFractalFriction::* ScalarFieldFunction)() const; 
  typedef Vec3 (CFractalFriction::* VectorFieldFunction)() const; 
	typedef pair<bool,double> (CFractalFriction::* CheckedScalarFieldFunction)() const;
	
  // type & dummy implementation for parameter setting function 
  typedef void (CFractalFriction::* ScalarSetFunction)(double);
  static ScalarSetFunction getScalarSetFunction(const string&){return NULL;};

 private:

 public:
  CFractalFriction();
  CFractalFriction(CParticle*,CParticle*,const FractalFrictionIGP&);
  virtual ~CFractalFriction();
  
  static string getType() {return "FractalFriction";};
  
  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);
  static CheckedScalarFieldFunction getCheckedScalarFieldFunction(const string&);
  friend ostream& operator<<(ostream&,const CFractalFriction&);
  friend class TML_PackedMessageInterface;
};

#endif //__FRACTALFRICTION_H
