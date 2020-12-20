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


#ifndef MODEL_LOCALDAMPING_HPP
#define MODEL_LOCALDAMPING_HPP

#include "Model/LocalDamping.h"

/*!
  Construct a local damping "interaction" for a particle

  \param P the particle
  \param visc the damping coefficient
  \param dt the time step
*/
template <class T>
CLocalDamping<T>::CLocalDamping(T* P,double visc, double dt)
{
  m_p=P;
  m_visc=visc;
  m_dt=dt;
  m_force=Vec3(0.0,0.0,0.0);
}

/*!
  Construct a local damping "interaction" for a particle

  \param P the particle
  \param param the parameters
*/
template <class T>
CLocalDamping<T>::CLocalDamping(T* P,const CLocalDampingIGP& param)
{
  m_p=P;
  m_visc=param.getVisc();
  m_dt=param.getTimeStep();
  m_force=Vec3(0.0,0.0,0.0);
}

/*!
  Construct a local damping "interaction" for a particle

  \param P the particle
  \param param pointer to the parameters
*/
template <class T>
CLocalDamping<T>::CLocalDamping(T* P,CLocalDampingIGP* param)
{
  m_p=P;
  m_visc=param->getVisc();
  m_dt=param->getTimeStep();
  m_force=Vec3(0.0,0.0,0.0);
}

/*!
  destructor
*/
template <class T>
CLocalDamping<T>::~CLocalDamping()
{}

template <class T>
void CLocalDamping<T>::setTimeStepSize(double dt)
{
  m_dt = dt;
}

/*!
  Calculate the local damping force.
*/
template <class T>
void CLocalDamping<T>::calcForces()
{
  m_E_diss=0.0;// zero dissipated energy
  Vec3 v=m_p->getVel();
  Vec3 frc=m_p->getForce();

  double dampFx, dampFy, dampFz;

  dampFx = fabs(frc.X());
  if (v.X() < 0) {
     dampFx *= -1.;
  }
  if (v.X() == 0) {
     dampFx = 0.0;
  }
  dampFy = fabs(frc.Y());
  if (v.Y() < 0) {
     dampFy *= -1.;
  }
  if (v.Y() == 0) {
     dampFy = 0.0;
  }
  dampFz = fabs(frc.Z());
  if (v.Z() < 0) {
     dampFz *= -1.;
  }
  if (v.Z() == 0) {
     dampFz = 0.0;
  }

  m_force=-1.0*m_visc*Vec3(dampFx, dampFy, dampFz);

  m_p->applyForce(m_force,m_p->getPos());  

  m_E_diss=m_visc*m_force.norm()*v.norm()*m_dt;
}

/*!
  Get the particle member function which returns a scalar field of a given name.
 
  \param name the name of the field
*/
template <class T>
typename CLocalDamping<T>::ScalarFieldFunction CLocalDamping<T>::getScalarFieldFunction(const string& name)
{
  typename CLocalDamping<T>::ScalarFieldFunction sf;

  if(name=="dissipated_energy"){
    sf=&CLocalDamping<T>::getDissipatedEnergy;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction scalar  access function" << endl;
  }

  return sf;
}

/*!
  Get the particle member function which returns a checked scalar field of a given name.

  \param name the name of the field
*/
template <class T>
typename CLocalDamping<T>::CheckedScalarFieldFunction CLocalDamping<T>::getCheckedScalarFieldFunction(const string& name)
{
  typename CLocalDamping<T>::CheckedScalarFieldFunction sf;

  sf=NULL;
  cerr << "ERROR - invalid name for interaction scalar  access function" << endl;
  
  return sf;
}


/*!
  Get the particle member function which returns a vector field of a given name.
 
  \param name the name of the field
*/
template <class T>
typename CLocalDamping<T>::VectorFieldFunction CLocalDamping<T>::getVectorFieldFunction(const string& name)
{
  typename CLocalDamping<T>::VectorFieldFunction vf;

  if (name=="force"){
    vf=&CLocalDamping<T>::getForce;
  } else {
      vf=NULL;
      cerr << "ERROR - invalid name for interaction vector access function" << endl;
  }
  
  return vf;
}

/*!
  return the amount of energy dissipated during the last time step
*/
template <class T>
double CLocalDamping<T>::getDissipatedEnergy() const
{
  return m_E_diss;
}

template <class T>
Vec3 CLocalDamping<T>::getForce() const
{
  return m_force;
}

/*!
  check if any of the particles in the interaction fits tag & mask

  \param tag the tag
  \param mask the mask
*/
template <class T>
bool CLocalDamping<T>::hasTag(int tag,int mask) const
{
 int tag1=m_p->getTag();

  return ((tag1 & mask)==(tag & mask));
}

/*!
  return a vector of all particle IDs
*/
template <class T>
vector<int> CLocalDamping<T>::getAllID() const
{
  vector<int> res;

  res.push_back(m_p->getID());

  return res;
}

#endif
