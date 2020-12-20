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


#ifndef MODEL_DAMPING_HPP
#define MODEL_DAMPING_HPP

#include "Model/Damping.h"
//#include "Foundation/console.h"

template <class T>
double CDamping<T>::s_limit2=1e-12; // default error limit : 1e-6 
template <class T>
int CDamping<T>::s_flops = 0;


/*!
  Construct a damping "interaction" for a particle

  \param P the particle
  \param V the reference velocity
  \param visc the artificial viscosity
  \param dt the time step 
  \param mi the maximum number of iterations
*/
template <class T>
CDamping<T>::CDamping(T* P,const Vec3& V,double visc,double dt,int mi)
{
  m_p=P;
  m_vref=V;
  m_visc=visc;
  m_dt=dt;
  m_maxiter=mi;
  m_force=Vec3(0.0,0.0,0.0);
}

/*!
  Construct a damping "interaction" for a particle

  \param P the particle
  \param param the parameters
*/
template <class T>
CDamping<T>::CDamping(T* P,const CDampingIGP& param)
{
  m_p=P;
  m_vref=param.getVRef();
  m_visc=param.getVisc();
  m_dt=param.getTimeStep();
  m_maxiter=param.getMaxIter();
  m_force=Vec3(0.0,0.0,0.0);
}

/*!
  Construct a damping "interaction" for a particle

  \param P the particle
  \param param pointer to the parameters
*/
template <class T>
CDamping<T>::CDamping(T* P,CDampingIGP* param)
{
  m_p=P;
  m_vref=param->getVRef();
  m_visc=param->getVisc();
  m_dt=param->getTimeStep();
  m_maxiter=param->getMaxIter();
  m_force=Vec3(0.0,0.0,0.0);
}

/*!
  destructor
*/
template <class T>
CDamping<T>::~CDamping()
{}

template <class T>
void CDamping<T>::setTimeStepSize(double dt)
{
  m_dt = dt;
}

/*!
  Calculate the damping force.

  25*count+8 flops
*/
template <class T>
void CDamping<T>::calcForces()
{
  m_E_diss=0.0;// zero dissipated energy
  Vec3 v=m_p->getVel();
  Vec3 v_rel=Vec3(0.0,0.0,0.0);
  Vec3 frc=m_p->getForce();
  
  double s=m_p->getInvMass();     
  double mass=m_p->getMass();
  double error=1.0;
  int count=0;
  while((error*error>s_limit2) & (count<m_maxiter)){ // 1 flop
    count++;
    Vec3 v_old=v_rel;
    v_rel=v-m_vref+s*m_dt*(frc-v_rel*m_visc*mass);        // 16 flops
    error=(v_rel-v_old).norm2();                     // 8 flops  
  }
  if(count>=m_maxiter){
    //console.Warning() << "damping diverges for particle  " << m_p->getID() << "error= " << error << "\n";
    v_rel=Vec3(0.0,0.0,0.0);
  }
  m_force=-1.0*m_visc*v_rel*mass;
  m_p->applyForce(m_force,m_p->getPos());  // 3 flops
  // m_E_diss=0.5*(v*v-(v_rel+m_vref)*(v_rel+m_vref))*m_p->getMass();
  m_E_diss=m_visc*v_rel*v*m_dt;
}

/*!
  Get the particle member function which returns a scalar field of a given name.
 
  \param name the name of the field
*/
template <class T>
typename CDamping<T>::ScalarFieldFunction CDamping<T>::getScalarFieldFunction(const string& name)
{
  typename CDamping<T>::ScalarFieldFunction sf;

  if(name=="dissipated_energy"){
    sf=&CDamping<T>::getDissipatedEnergy;
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
typename CDamping<T>::CheckedScalarFieldFunction CDamping<T>::getCheckedScalarFieldFunction(const string& name)
{
  typename CDamping<T>::CheckedScalarFieldFunction sf;

  sf=NULL;
  cerr << "ERROR - invalid name for interaction scalar  access function" << endl;
  
  return sf;
}


/*!
  Get the particle member function which returns a vector field of a given name.
 
  \param name the name of the field
*/
template <class T>
typename CDamping<T>::VectorFieldFunction CDamping<T>::getVectorFieldFunction(const string& name)
{
  typename CDamping<T>::VectorFieldFunction vf;

  if (name=="force"){
    vf=&CDamping<T>::getForce;
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
double CDamping<T>::getDissipatedEnergy() const
{
  return m_E_diss;
}

template <class T>
Vec3 CDamping<T>::getForce() const
{
  return m_force;
}

/*!
  check if any of the particles in the interaction fits tag & mask

  \param tag the tag
  \param mask the mask
*/
template <class T>
bool CDamping<T>::hasTag(int tag,int mask) const
{
 int tag1=m_p->getTag();

  return ((tag1 & mask)==(tag & mask));
}

/*!
  return a vector of all particle IDs
*/
template <class T>
vector<int> CDamping<T>::getAllID() const
{
  vector<int> res;

  res.push_back(m_p->getID());

  return res;
}

#endif
