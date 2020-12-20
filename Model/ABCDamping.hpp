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

#ifndef MODEL_ABCDAMPING_HPP
#define MODEL_ABCDAMPING_HPP

#include<cmath>
using std::fabs;
using std::exp;



/*!
 Construct a damping "interaction" for a particle

  \param P the particle
  \param param the parameters
*/
template <class ParticleType>
ABCDamping<ParticleType>::ABCDamping(ParticleType* P, ABCDampingIGP* param) 
: CDamping<ParticleType>(P,param)
{
  // calc local viscosity
  Vec3 pos=param->getPos();
  Vec3 normal=param->getNormal();
  double c1=param->getC1();
  double v_0=this->m_visc;
  double dist=fabs((this->m_p->getInitPos()-pos)*normal);
  double v_new=v_0/exp(dist*c1);
  this->m_visc=v_new;
}

template <class ParticleType>
ABCDamping<ParticleType>::~ABCDamping()
{}


/*  Get the particle member function which returns a scalar field of a given name.
 
  \param name the name of the field
*/
template <class T>
typename ABCDamping<T>::ScalarFieldFunction ABCDamping<T>::getScalarFieldFunction(const string& name)
{
  typename ABCDamping<T>::ScalarFieldFunction sf;
  
  sf=NULL;
  cerr << "ERROR - invalid name for interaction scalar  access function" << endl;
  
  return sf;
}


/*!
  Get the particle member function which returns a checked scalar field of a given name.

  \param name the name of the field
*/
template <class T>
typename ABCDamping<T>::CheckedScalarFieldFunction ABCDamping<T>::getCheckedScalarFieldFunction(const string& name)
{
  typename ABCDamping<T>::CheckedScalarFieldFunction sf;

  sf=NULL;
  cerr << "ERROR - invalid name for interaction scalar  access function" << endl;
  
  return sf;
}


/*!
  Get the particle member function which returns a vector field of a given name.
 
  \param name the name of the field
*/
template <class T>
typename ABCDamping<T>::VectorFieldFunction ABCDamping<T>::getVectorFieldFunction(const string& name)
{
  typename ABCDamping<T>::VectorFieldFunction vf;

  vf=NULL;
  cerr << "ERROR - invalid name for interaction vector access function" << endl;
  
  return vf;
}

#endif //MODEL_ABCDAMPING_HPP
