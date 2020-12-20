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

#ifndef __ROT_DAMPING_HPP
#define __ROT_DAMPING_HPP

template <class T>
double CRotDamping<T>::s_limit2=1e-12; // default error limit : 1e-6 
template <class T>
int CRotDamping<T>::s_flops = 0;


/*!
  Construct a rotational damping "interaction" for a particle

  \param P the particle
  \param param pointer to the parameters
*/
template <class T>
CRotDamping<T>::CRotDamping(T* P,CDampingIGP* param)
{
  m_p=P;
  m_vref=param->getVRef();
  m_visc=param->getVisc();
  m_dt=param->getTimeStep();
  m_maxiter=param->getMaxIter();
}

/*!
  destructor
*/
template <class T>
CRotDamping<T>::~CRotDamping()
{
}

template <class T>
void CRotDamping<T>::setTimeStepSize(double dt)
{
  m_dt = dt;
}

/*!
  Calculate the damping force.

  25*count+8 flops
*/
template <class T>
void CRotDamping<T>::calcForces()
{
  m_E_diss=0.0;// zero dissipated energy
  Vec3 v=m_p->getAngVel();
  Vec3 v_rel=Vec3(0.0,0.0,0.0);
  Vec3 frc=m_p->getMoment();
  
  double s=1.0/m_p->getInertRot();       
  double in=m_p->getInertRot();
  double mass=m_p->getMass();

  double error=1.0;
  int count=0;
  while((error*error>s_limit2) & (count<m_maxiter)){ // 1 flop
    count++;
    Vec3 v_old=v_rel;
    v_rel=v-m_vref+s*m_dt*(frc-v_rel*m_visc*in);        // 16 flops
    error=(v_rel-v_old).norm2();                     // 8 flops  
  }
  if(count>=m_maxiter){
    //console.Warning() << "damping diverges for particle  " << m_p->getID() << "error= " << error << "\n";
    v_rel=Vec3(0.0,0.0,0.0);
  }
  m_force=-1.0*m_visc*v_rel*mass;
  m_p->applyMoment(-1.0*m_visc*v_rel*in);  // 3 flops
  m_E_diss=m_visc*v_rel*v*m_dt; // wrong
}

/*!
  Get the particle member function which returns a scalar field of a given name.
 
  \param name the name of the field
*/
template <class T>
typename CRotDamping<T>::ScalarFieldFunction CRotDamping<T>::getScalarFieldFunction(const string& name)
{
  typename CRotDamping<T>::ScalarFieldFunction sf;

  if(name=="dissipated_energy"){
    sf=&CRotDamping<T>::getDissipatedEnergy;
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
typename CRotDamping<T>::CheckedScalarFieldFunction CRotDamping<T>::getCheckedScalarFieldFunction(const string& name)
{
  typename CRotDamping<T>::CheckedScalarFieldFunction sf;

  sf=NULL;
  cerr << "ERROR - invalid name for interaction scalar  access function" << endl;
  
  return sf;
}


/*!
  Get the particle member function which returns a vector field of a given name.
 
  \param name the name of the field
*/
template <class T>
typename CRotDamping<T>::VectorFieldFunction CRotDamping<T>::getVectorFieldFunction(const string& name)
{
  typename CRotDamping<T>::VectorFieldFunction vf;

  if (name=="force"){
    vf=&CRotDamping<T>::getForce;
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
double CRotDamping<T>::getDissipatedEnergy() const
{
  return m_E_diss;
}

template <class T>
Vec3 CRotDamping<T>::getForce() const
{
  return m_force;
}

/*!
  check if any of the particles in the interaction fits tag & mask

  \param tag the tag
  \param mask the mask
*/
template <class T>
bool CRotDamping<T>::hasTag(int tag,int mask) const
{
 int tag1=m_p->getTag();

  return ((tag1 & mask)==(tag & mask));
}

/*!
  return a vector of all particle IDs
*/
template <class T>
vector<int> CRotDamping<T>::getAllID() const
{
  vector<int> res;

  res.push_back(m_p->getID());

  return res;
}

#endif // __ROT_DAMPING_HPP
