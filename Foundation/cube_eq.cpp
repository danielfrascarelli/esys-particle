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

#include "cube_eq.h"

// --- system includes ---
#include <cmath>
using std::sqrt;

/*!
  construct cubic equation of the form x^3+ax^2+bx+c

  \param a coefficient for x^2
  \param b coefficient for x
  \param c constant coefficient
*/
CubicEquation::CubicEquation(double a,double b,double c)
{
  m_a=a;
  m_b=b;
  m_c=c;
}

/*!
  Get the roots. Get one root (r_1) by a bisection method and the other 2 (if real) by solving the 
  quadratic equation resulting from dividing the eqation by (x-r_1). Returns the roots as a STL-set 
  so they are ordered.

  \param tol the precision of the calculation
  \param valid returns the validity of the result, i.e. if valid==false there was no positive root found
*/
set<double> CubicEquation::getRealRoots(double tol)
{
  set<double> rootset;
  double root1;

  double x_0=-1.0*sqrt(1+m_a*m_a+m_b*m_b+m_c*m_c);
  double x_1=-1.0*x_0;
  
  if((f(x_0)*f(x_1))<0.0){ // change of sign within range
    // get one root
    root1=bisect(x_0,x_1,tol);
    rootset.insert(root1);
    // divide eq by (x-root1) => quadratic eq x^2+px+q
    double p=m_a+root1;
    double q=m_b+(m_a+root1)*root1;
    // solve quadratic eq.
    double d=0.25*p*p-q; 
    if(d>=0.0){ // 2 real roots
      double dsr=sqrt(d);
      rootset.insert(-0.5*p-dsr);
      rootset.insert(-0.5*p+dsr);
    }
  } else {
    // throw some error
    // can't happen -> cubic eq has always at least 1 real root
  }

  return rootset;
}

/*!
  return x^3+ax^2+bx+c

  \param x
*/
double CubicEquation::f(double x) 
{
  return x*x*x+m_a*x*x+m_b*x+m_c;
}

double CubicEquation::bisect(double x_0,double x_1,double tol)
{
  double res;

  if(x_1-x_0<tol){
    res=0.5*(x_0+x_1);
  } else{
    double x_h=0.5*(x_0+x_1);
    if(f(x_0)*f(x_h)<0.0){ // change of sign in lower half
      res=bisect(x_0,x_h,tol);
    } else { // change of sign must be in upper half
      res=bisect(x_h,x_1,tol);
    }
  }
  return res;
}
