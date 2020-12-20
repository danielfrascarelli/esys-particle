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

#include <mpi.h>
#include "Model/RotThermParticle.h"
#include "Geometry/SimpleParticleData.h"
#include "Foundation/console.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

CRotThermParticle::CRotThermParticle() : CRotParticleVi(),CThermParticle()
{
 // m_global_id=-1;
//  flag=false;
}


CRotThermParticle::CRotThermParticle(
  const esys::lsm::SimpleParticleData &particleData
)
  : CRotParticleVi(particleData),
    CThermParticle(particleData.getRadius())
{// ????
}

CRotThermParticle::CRotThermParticle(const CRotParticleVi &p)
  : CRotParticleVi(p),
    CThermParticle(p.getRad())
{
}

CRotThermParticle::CRotThermParticle(const CParticle &p)
  : CRotParticleVi(p),
    CThermParticle(p.getRad())
{
}

CRotThermParticle::CRotThermParticle(
  double rad,
  double mass,
  const Vec3& pos,
  const Vec3& vel,
  const Vec3& force,
  int id,
  bool is_dyn
) :
    CRotParticleVi(rad, mass, pos, vel, force, id, is_dyn),
    CThermParticle(rad)
{
}

/*!
  Construct particle. Old and initial position are assumed to be identical to current position. 

  \param rad radius
  \param mass mass
  \param pos current position
  \param vel current velocity
  \param force currently applied force
  \param id particle id
  \param quat particel quaternion
  \param inertRot inert of rotation
  \param moment currently applied tarque
  \param angvel current angular velocity  
*/

CRotThermParticle::CRotThermParticle(
  double rad,
  double mass,
  const Vec3& pos,
  const Vec3& vel,
  const Vec3& force,
  int id,
  Quaternion & quat,
  double inertRot,
  const Vec3& moment,
  const Vec3& angvel,
  const Vec3& angvel_t,
  double temperature,
  double temperature_ini,
  double Cp,
  double heat_frict,
  double heat_trans,
  double therm_expansion0,
  double therm_expansion1,
  double therm_expansion2
) :
  CRotParticleVi(
    rad,
    mass,
    pos,
    vel,
    force,
    id,
    quat,
    inertRot,
    moment,
    angvel,
    angvel_t
  ),
  CThermParticle(
    temperature,
    temperature_ini,
    Cp,
    heat_frict,
    heat_trans,
    therm_expansion0,
    therm_expansion1,
    therm_expansion2,
    rad
  )
{
}

CRotThermParticle::CRotThermParticle(
  double rad,
  double mass,
  const Vec3& pos,
  const Vec3& oldpos,
  const Vec3& initpos,
  const Vec3& vel,
  const Vec3& force,
  int id,
  const Quaternion & quat,
  const Quaternion & initquat,
  double inertRot,
  const Vec3& moment,
  const Vec3& angvel,
  const Vec3& angvel_t,
  double temperature,
  double temperature_ini,
  double Cp,
  double heat_frict,
  double heat_trans,
  double therm_expansion0,
  double therm_expansion1,
  double therm_expansion2
)
  :
    CRotParticleVi(
      rad,
      mass,
      pos,
      oldpos,
      initpos,
      vel,
      force,
      id,
      quat,
      initquat,
      inertRot,
      moment,
      angvel,
      angvel_t
  ),
  CThermParticle(
    temperature,
    temperature_ini,
    Cp,
    heat_frict,
    heat_trans,
    therm_expansion0,
    therm_expansion1,
    therm_expansion2,
    rad
  )
{
}


/*!
  Pack a CParticle into a TML packed message

  \param p the particle
  \todo BasicParticle data should be handled by pack<Basicparticle>
*/
template<>
void TML_PackedMessageInterface::pack<CRotThermParticle>(const CRotThermParticle& p)
{
  append(p.m_tag);
  append(p.m_pos.X()) ;
  append(p.m_pos.Y()) ;
  append(p.m_pos.Z()) ;
  append(p.m_oldpos.X());
  append(p.m_oldpos.Y());
  append(p.m_oldpos.Z());
  append(p.m_initpos.X());
  append(p.m_initpos.Y());
  append(p.m_initpos.Z());
  append(p.m_circular_shift.X());
  append(p.m_circular_shift.Y());
  append(p.m_circular_shift.Z());
  append(p.m_vel.X());
  append(p.m_vel.Y());
  append(p.m_vel.Z());
  append(p.m_force.X());
  append(p.m_force.Y());
  append(p.m_force.Z());
  append(p.m_rad);
  append(p.m_mass);
//  append(p.m_global_id);
 
// wyc added !

  append(p.m_inertRot);
  append(p.m_moment.X()) ;
  append(p.m_moment.Y()) ;
  append(p.m_moment.Z()) ;
  append(p.m_angVel.X()) ;
  append(p.m_angVel.Y()) ;
  append(p.m_angVel.Z()) ;  
  append(p.m_quat.return_sca()) ;
  append(p.m_quat.return_vec().X()) ;
  append(p.m_quat.return_vec().Y()) ;
  append(p.m_quat.return_vec().Z()) ; 
  append(p.m_initquat.return_sca()) ;
  append(p.m_initquat.return_vec().X()) ;
  append(p.m_initquat.return_vec().Y()) ;
  append(p.m_initquat.return_vec().Z()) ; 
  append(p.m_angVel_t.X()) ;
  append(p.m_angVel_t.Y()) ;
  append(p.m_angVel_t.Z()) ;

  append(p.m_temperature);
  append(p.m_temperature_ini);
  append(p.m_Cp) ;
  append(p.m_heat_frict);
  append(p.m_heat_trans); 
  append(p.m_therm_expansion0);
  append(p.m_therm_expansion1);
  append(p.m_therm_expansion2);
  append(p.m_rad_ini);

  append(p.m_global_id); // original one here
}


/*!
  Unpack a CParticle from a TML packed message

  \param p the particle
*/

template<>
void TML_PackedMessageInterface::unpack<CRotThermParticle>(CRotThermParticle& p)
{
  const int numElems = 47 ; //????
  double db[numElems] ;

  p.m_tag=pop_int();
  pop_doubles(db, numElems);
  p.m_pos            = Vec3(db[0],db[1],db[2]) ;
  p.m_oldpos         = Vec3(db[3],db[4],db[5]) ;
  p.m_initpos        = Vec3(db[6],db[7],db[8]) ;
  p.m_circular_shift = Vec3(db[9],db[10],db[11]) ;
  p.m_vel            = Vec3(db[12],db[13],db[14]) ;
  p.m_force          = Vec3(db[15],db[16],db[17]) ;
  p.m_rad            = db[18];
  p.m_mass           = db[19];
  p.m_div_mass = 1.0/p.m_mass;

  p.m_inertRot   = db[20];
  p.m_div_inertRot = (p.m_inertRot!=0.0) ? 1.0/p.m_inertRot : 0.0;
  p.m_moment     = Vec3(db[21],db[22],db[23]);
  p.m_angVel     = Vec3(db[24],db[25],db[26]);
  p.m_quat       = Quaternion(db[27],Vec3(db[28],db[29],db[30])); 
  p.m_initquat   = Quaternion(db[31],Vec3(db[32],db[33],db[34]));
  p.m_angVel_t = Vec3(db[35],db[36],db[37]) ;

  p.m_temperature = db[38] ;
  p.m_temperature_ini = db[39] ;
  p.m_Cp = db[40];
  p.m_heat_frict = db[41];
  p.m_heat_trans= db[42];
  p.m_therm_expansion0 =db[43] ;
  p.m_therm_expansion1 =db[44] ;
  p.m_therm_expansion2 =db[45] ;
  p.m_rad_ini = db[46] ;
 
  p.m_global_id=pop_int();
}

  
map<string,AField*> CRotThermParticle::generateFields(ParallelParticleArray<CRotThermParticle>* particles)
{
  map<string,AField*> res;

  return res;
}



CRotThermParticle::exchangeType CRotThermParticle::getExchangeValues()
{
  return 
    exchangeType(
      m_pos     - m_circular_shift,
      m_initpos - m_circular_shift,
      m_vel,
      m_angVel,
      m_angVel_t,
      m_quat,
      m_temperature,
      m_temperature_ini
    );
}

/*!
  Set pos, vel and angular vel from exchangeType
                                                                                
  \param E the exchanged values
*/
void CRotThermParticle::setExchangeValues(const exchangeType& E)
{
   m_pos        = E.m_pos     + m_circular_shift;
   m_initpos    = E.m_initPos + m_circular_shift;
   m_vel        = E.m_vel;
   m_angVel = E.m_angVel;
   m_angVel_t = E.m_angVel_t;
   m_quat       = E.m_quat;
   m_temperature= E.m_temperature;
   m_temperature_ini= E.m_temperature_ini;
}

template<>
void TML_PackedMessageInterface::pack<CRotThermParticle::exchangeType>(const CRotThermParticle::exchangeType& p)
{
  append(p.m_pos.X());
  append(p.m_pos.Y());
  append(p.m_pos.Z());
  append(p.m_initPos.X());
  append(p.m_initPos.Y());
  append(p.m_initPos.Z());
  append(p.m_vel.X());
  append(p.m_vel.Y());
  append(p.m_vel.Z());
  append(p.m_angVel.X());
  append(p.m_angVel.Y());
  append(p.m_angVel.Z());
  append(p.m_angVel_t.X());
  append(p.m_angVel_t.Y());
  append(p.m_angVel_t.Z());  
  append(p.m_quat.return_sca());
  append(p.m_quat.return_vec().X());
  append(p.m_quat.return_vec().Y());
  append(p.m_quat.return_vec().Z());
  append(p.m_temperature) ;
  append(p.m_temperature_ini) ;
}

/*!
  Unpack an exchangeType from a TML packed message
*/
template<>
void TML_PackedMessageInterface::unpack<CRotThermParticle::exchangeType>(CRotThermParticle::exchangeType& p)
{
  const int numDoubles = 21 ;
  double db[numDoubles];

  pop_doubles(db, numDoubles);
  p.m_pos        = Vec3(db[0],db[1],db[2]);
  p.m_initPos    = Vec3(db[3],db[4],db[5]);
  p.m_vel        = Vec3(db[6],db[7],db[8]);
  p.m_angVel     = Vec3(db[9],db[10],db[11]);
  p.m_angVel_t     = Vec3(db[12],db[13],db[14]);
  p.m_quat       = Quaternion(db[15],Vec3(db[16],db[17],db[18]));
  p.m_temperature = db[19] ;
  p.m_temperature_ini = db[20] ;
}

CRotThermParticle::ScalarFieldFunction CRotThermParticle::getScalarFieldFunction(const string& name)
{
  CRotThermParticle::ScalarFieldFunction sf;
                                                                                
  if(name=="id"){
    sf=&CRotThermParticle::getIDField;
  } else if(name=="sigma_xx_2d"){
    sf=&CRotThermParticle::sigma_xx_2D;
  } else if(name=="sigma_xy_2d"){
    sf=&CRotThermParticle::sigma_xy_2D;
  } else if(name=="sigma_yy_2d"){
    sf=&CRotThermParticle::sigma_yy_2D;
  } else if(name=="sigma_d"){
    sf=&CRotThermParticle::sigma_d;
  } else if(name=="sigma_vonMises"){
    sf=&CRotThermParticle::sigma_vonMises;
  } else if(name=="e_kin"){
    sf=&CRotThermParticle::getKineticEnergy;
  } else if(name=="radius"){
    sf=&CRotThermParticle::getRad;
  } else if(name=="temperature"){
    sf=&CRotThermParticle::getTemperature;

  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for particle scalar  access function" << endl;                                                                                
  }

  return sf;
}


CRotThermParticle::VectorFieldFunction CRotThermParticle::getVectorFieldFunction(const string& name)
{
  CRotThermParticle::VectorFieldFunction sf;
                                                                                
  if(name=="displacement"){
    sf=&CRotThermParticle::getTotalDisplacement;
  } else if(name=="velocity"){
    sf=&CRotThermParticle::getVel;
  } else if(name=="ang_velocity"){
    sf=&CRotThermParticle::getAngVelNR;
  } else if(name=="position"){
    sf=&CRotThermParticle::getPos;
  } else if(name=="force"){
    sf=&CRotThermParticle::getForce;
  } else if(name=="anglevector"){
    sf=&CRotParticleVi::getAngVector;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for particle vector access function" << endl;
  }
                                                                                
  return sf;
}

// used in 2D to output rotated angle
/*
Vec3 CRotParticleVi::getAngVector() const
{
 double half_cos, half_sin,s_in, c_os  ;
 Vec3 AngVector ;
  half_sin = getQuat().return_vec().Z() ;
  half_cos = getQuat().return_sca();
 s_in = 2.0*half_cos*half_sin ;
 c_os = half_cos*half_cos -half_sin*half_sin ;
 AngVector = Vec3(getRad()*c_os,getRad()*s_in,0.0) ;

 return AngVector ;
}
*/

//3D case: direction--- q1,q2,q3, size --- angle
/*
Vec3 CRotParticleVi::getAngVector() const
{
 double half_cos, half_sin, c_os, angle  ;
 Vec3 AngVector,direction ;
  half_cos = getQuat().return_sca();
  half_sin = sqrt(abs(1-half_cos*half_cos));
  c_os = half_cos*half_cos -half_sin*half_sin ;
  if(c_os>=1.0) angle =0.0;
  else angle = acos(c_os);

 direction = Vec3(getQuat().return_vec().X(),getQuat().return_vec().Y(),getQuat().return_vec().Z() ) ;

 if(direction.norm()==0.0) direction = Vec3(0.0,0.0,0.0) ;
 else  direction = direction/direction.norm() ;

  AngVector = angle*direction ;

 return AngVector ;
}
*/




/*!
  set circular shift vector
                                                                                
  \param cv the circular shift vector
*/
void CRotThermParticle::setCircular(const Vec3& cv)
{
  CRotParticleVi::setCircular(cv);
}

/**
 * Save check-point data.
 */
/*void CRotThermParticle::saveCheckPointData(std::ostream& oStream)
{
  CParticle::saveCheckPointData(oStream);
  const char delim = ' ';
  oStream
    << delim
    << getQuat()
    << delim
    << getAngVel();
}

void CRotThermParticle::loadCheckPointData(std::istream &iStream)
{
  CParticle::loadCheckPointData(iStream);
  iStream >> m_quat;
  iStream >> m_angVel_t;
}
*/

void CRotThermParticle::zeroHeat()
{
  m_heat_frict = 0.0 ;
  m_heat_trans = 0.0 ;

}

void CRotThermParticle::applyHeatTrans(const double heat_trans)
{
  m_heat_trans += heat_trans ; 

}

void CRotThermParticle::applyHeatFrict(const double heat_frict)
{
  m_heat_frict += heat_frict ;

}



void CRotThermParticle::integrateTherm(double dt)
{
  m_temperature += m_heat_frict/(getMass()*getCp()) + dt*m_heat_trans ;
}

void CRotThermParticle::integrate(double dt)
{
  CRotParticleVi::integrate(dt) ;
}

void CRotThermParticle::thermExpansion() 
{
#if 0
  /*
   * non-linear expansion.
   */
  m_rad =
    m_rad_ini*(
      1.0 + 
      m_therm_expansion0 +
      m_therm_expansion1*(m_temperature- m_temperature_ini) +
      m_therm_expansion2*
        (m_temperature- m_temperature_ini)*
        (m_temperature- m_temperature_ini)
    ) ;
#else
  /*
   * linear expansion.
   */

  m_rad =
    m_rad_ini*
    (
      1.0
      +
      m_therm_expansion1*(m_temperature- m_temperature_ini)
    );
#endif
}

ostream& operator<<(ostream& ost, const CRotThermParticle& CP)
{//????
  ost << "--CParticle " << CP.m_global_id << "  --\n";
  ost << "Radius : " << CP.m_rad << " Mass : " << CP.m_mass << "\n";
  ost << "Position : " << CP.m_pos << "\n";
  ost << "Velocity : " << CP.m_vel << "\n";
  ost << "Force    : " << CP.m_force << "\n";
  return ost;
}

