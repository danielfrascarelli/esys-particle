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


#include "Model/RotParticle.h"
#include "Geometry/SimpleParticleData.h"
#include "Foundation/console.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

CRotParticle::CRotParticle() : CParticle()
{
  m_global_id=-1;
  flag=false;
}

CRotParticle::CRotParticle(const esys::lsm::SimpleParticleData &particleData) : CParticle(particleData)
{
  m_quat = Quaternion(1.0,Vec3(0.0,0.0,0.0));
  m_initquat = m_quat;
  if(getDo2dCalculations()){
    m_inertRot  = 0.5*particleData.getMass()*particleData.getRadius()*particleData.getRadius();// cylinder
  } else {
    m_inertRot  = 0.4*particleData.getMass()*particleData.getRadius()*particleData.getRadius();// sphere
  }
  m_div_inertRot = 1.0/m_inertRot;
  m_angVel = Vec3(0.0,0.0,0.0);
  m_moment     = Vec3(0.0,0.0,0.0);
  m_is_rot=true; // particleData has no is_rot info -> default to rot. dynamics on
}

CRotParticle::CRotParticle(const CParticle &particle) : CParticle(particle)
{
  m_quat = Quaternion(1.0,Vec3(0.0,0.0,0.0));
  m_initquat = m_quat;
  if(getDo2dCalculations()){
    m_inertRot  = 0.5*particle.getMass()*particle.getRad()*particle.getRad(); // cylinder
  } else {
    m_inertRot  = 0.4*particle.getMass()*particle.getRad()*particle.getRad(); // sphere
  }
  m_div_inertRot = 1.0/m_inertRot;
  m_angVel = Vec3(0.0,0.0,0.0);
  m_moment     = Vec3(0.0,0.0,0.0);
  m_is_rot=true; // CParticle has no is_rot info -> default to rot. dynamics on
}

/*!
  Construct particle with default rotational orientation. Inertia 
  is calculated from mass & radius. Used from Python interface.
*/
CRotParticle::CRotParticle(
  double rad,
  double mass,
  const Vec3& pos,
  const Vec3& vel,
  const Vec3& force,
  int id,
  bool is_dyn,
  bool is_rot
  ) : CParticle(rad, mass, pos, vel, force, id, is_dyn)
{
  m_quat = Quaternion(1.0,Vec3(0.0,0.0,0.0));
  m_initquat = m_quat;
  if(getDo2dCalculations()){
    m_inertRot  = 0.5*mass*rad*rad; // 2D -> cylinder
  } else {
    m_inertRot  = 0.4*mass*rad*rad; // 3D -> sphere
  }
  m_div_inertRot = 1.0/m_inertRot;
  m_angVel = Vec3(0.0,0.0,0.0);
  m_moment     = Vec3(0.0,0.0,0.0);
  m_is_rot=is_rot;
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
  \param is_rot rotational dynamics on/off
*/
CRotParticle::CRotParticle(
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
  bool is_dyn,
  bool is_rot
  ) :CParticle(rad, mass,pos,vel,force,id,is_dyn)
{
  m_circular_shift=Vec3(0.0,0.0,0.0);
  flag=false;

  m_quat         = quat;
  m_initquat     = quat;
  m_inertRot     = inertRot;
  m_moment       = moment;
  m_angVel       = angvel;
  m_div_inertRot = 1.0/m_inertRot;
  m_is_rot=is_rot;
}

CRotParticle::CRotParticle(
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
  bool is_dyn,
  bool is_rot
  ) :CParticle(rad, mass,pos,oldpos, initpos,vel,force,id,is_dyn)
{
  m_circular_shift=Vec3(0.0,0.0,0.0);
  flag=false;

  m_quat         = quat;
  m_initquat     = initquat;
  m_inertRot     = inertRot ;
  m_moment       = moment;
  m_angVel       = angvel;
  m_div_inertRot = 1.0/m_inertRot ;
  m_is_rot=is_rot;
}

/*!
  set the density of the particle
*/
void CRotParticle::setDensity(double rho)
{
  const double pi=3.141592654;
  
  if(getDo2dCalculations()){
    m_mass=rho*pi*m_rad*m_rad;
    m_inertRot  = 0.5*m_mass*m_rad*m_rad; // 2D -> cylinder
  } else {
    m_mass=(4.0/3.0)*rho*pi*m_rad*m_rad*m_rad;
    m_inertRot  = 0.4*m_mass*m_rad*m_rad; // 3D -> sphere	  
  }
  m_div_inertRot = 1.0/m_inertRot;
  if(m_mass!=0.0){
    m_div_mass=1.0/m_mass;
  } else {
    m_div_mass=0.0;
  }
}


/*!
  Pack a CParticle into a TML packed message

  \param p the particle
  \todo BasicParticle data should be handled by pack<Basicparticle>
*/
template<>
void TML_PackedMessageInterface::pack<CRotParticle>(const CRotParticle& p)
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
  
  append(p.m_global_id); // original one here
  
  append(p.m_is_dynamic);
  append(p.m_is_rot);
  append(p.m_div_vol);
}


/*!
  Unpack a CParticle from a TML packed message

  \param p the particle
*/

template<>
void TML_PackedMessageInterface::unpack<CRotParticle>(CRotParticle& p)
{
  const int numElems = 35;
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

  p.m_global_id=pop_int();
  p.m_is_dynamic=pop_bool();
  p.m_is_rot=pop_bool();
  p.m_div_vol=pop_double();
}

/*!
  Do the time integration for the particle.

  \param dt the time step
*/
void CRotParticle::integrate(double dt)
{
  if(m_is_rot){
    if (CParticle::getDo2dCalculations()) {
      m_force  = Vec3(m_force.X(), m_force.Y(), 0);
      m_moment = Vec3(0, 0, m_moment.Z());
    }

    // integrate rotational part
    m_angVel += (dt*m_div_inertRot) * m_moment;
    m_quat   += (dt/2.0)*(Quaternion(0, m_angVel)*m_quat);
  } else {
    m_angVel = Vec3(0.0,0.0,0.0);
  }
  // Integrate linear part
  CParticle::integrate(dt);
}

/*!
  zero forces on particle
*/
void CRotParticle::zeroForce()
{
  m_force = Vec3(0.0,0.0,0.0);
  m_moment = Vec3(0.0,0.0,0.0);
  m_sigma = Matrix3();
}

void CRotParticle::rescale()
{
  const double module = 
    sqrt(
      m_quat.return_sca()*m_quat.return_sca()
      + m_quat.return_vec().X()*m_quat.return_vec().X()
      + m_quat.return_vec().Y()*m_quat.return_vec().Y()
      + m_quat.return_vec().Z()*m_quat.return_vec().Z()
    );
    
  if(module != 0.0 ) {
    const double inverse = 1.0/module;
    m_quat  =  m_quat*inverse ;
  } else  {
    cerr << " Quaternion wrong !!!  " ;
  }
}

/*!
  Apply a moment to a particle at a given position. 

  \param moment  tarque  
*/

void CRotParticle::applyMoment(const Vec3& moment)
{
  m_moment += moment;
}

CRotParticle::exchangeType CRotParticle::getExchangeValues()
{
  return 
    exchangeType(
      m_pos     - m_circular_shift,
      m_initpos - m_circular_shift,
      m_vel,
      m_angVel,
      m_quat,
      m_is_dynamic,
      m_is_rot,
      m_sigma	
    );
}

/*!
  Set pos, vel and angular vel from exchangeType
                                                                                
  \param E the exchanged values
*/
void CRotParticle::setExchangeValues(const exchangeType& E)
{
  m_pos        = E.m_pos     + m_circular_shift;
  m_initpos    = E.m_initPos + m_circular_shift;
  m_vel        = E.m_vel;
  m_angVel = E.m_angVel;
  m_quat       = E.m_quat;
  m_is_dynamic = E.m_is_dynamic;
  m_is_rot = E.m_is_rot;
  m_sigma = E.m_stress;
}

template<>
void TML_PackedMessageInterface::pack<CRotParticle::exchangeType>(const CRotParticle::exchangeType& p)
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
  append(p.m_quat.return_sca());
  append(p.m_quat.return_vec().X());
  append(p.m_quat.return_vec().Y());
  append(p.m_quat.return_vec().Z());
  append(p.m_is_dynamic);
  append(p.m_is_rot);
  append(p.m_stress(0,0));
  append(p.m_stress(0,1));
  append(p.m_stress(0,2));
  append(p.m_stress(1,0));
  append(p.m_stress(1,1));
  append(p.m_stress(1,2));
  append(p.m_stress(2,0));
  append(p.m_stress(2,1));
  append(p.m_stress(2,2));
}

/*!
  Unpack an exchangeType from a TML packed message
*/
template<>
void TML_PackedMessageInterface::unpack<CRotParticle::exchangeType>(CRotParticle::exchangeType& p)
{
  const int numDoubles = 16;
  double db[numDoubles];

  pop_doubles(db, numDoubles);
  p.m_pos        = Vec3(db[0],db[1],db[2]);
  p.m_initPos    = Vec3(db[3],db[4],db[5]);
  p.m_vel        = Vec3(db[6],db[7],db[8]);
  p.m_angVel     = Vec3(db[9],db[10],db[11]);
  p.m_quat       = Quaternion(db[12],Vec3(db[13],db[14],db[15]));
  p.m_is_dynamic = pop_bool();
  p.m_is_rot     = pop_bool();
  
  double db2[9];
  pop_doubles(db2, 9); // stress tensor
  p.m_stress=Matrix3(db2);
  
}

CRotParticle::ScalarFieldFunction CRotParticle::getScalarFieldFunction(const string& name)
{
  CRotParticle::ScalarFieldFunction sf;
                                                                                
  if(name=="id"){
    sf=&CRotParticle::getIDField;
  } else if(name=="tag"){
    sf=&CRotParticle::getTagField;    
  } else if(name=="sigma_xx_2d"){
    sf=&CRotParticle::sigma_xx_2D;
  } else if(name=="sigma_xy_2d"){
    sf=&CRotParticle::sigma_xy_2D;
  } else if(name=="sigma_yy_2d"){
    sf=&CRotParticle::sigma_yy_2D;
  } else if(name=="sigma_d"){
    sf=&CRotParticle::sigma_d;
  } else if(name=="sigma_vonMises"){
    sf=&CRotParticle::sigma_vonMises;
  } else if(name=="e_kin"){
    sf=&CRotParticle::getKineticEnergy;
  } else if(name=="e_kin_rot"){
    sf=&CRotParticle::getAngularKineticEnergy;
  } else if(name=="e_kin_linear"){
    sf=&CRotParticle::getLinearKineticEnergy;
  } else if(name=="radius"){
    sf=&CRotParticle::getRad;
  } else if(name=="v_abs"){
    sf=&CRotParticle::getAbsVel;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for particle scalar  access function" << endl;                                                                                
  }

  return sf;
}


CRotParticle::VectorFieldFunction CRotParticle::getVectorFieldFunction(const string& name)
{
  CRotParticle::VectorFieldFunction sf;
                                                                                
  if(name=="displacement"){
    sf=&CRotParticle::getTotalDisplacement;
  } else if(name=="velocity"){
    sf=&CRotParticle::getVel;
  } else if(name=="ang_velocity"){
    sf=&CRotParticle::getAngVelNR;
  } else if(name=="position"){
    sf=&CRotParticle::getPos;
  } else if(name=="force"){
    sf=&CRotParticle::getForce;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for particle vector access function" << endl;
  }
                                                                                
  return sf;
}

/*!
  set circular shift vector
                                                                                
  \param cv the circular shift vector
*/
void CRotParticle::setCircular(const Vec3& cv)
{
  CParticle::setCircular(cv);
}

/*!
  Save snapshot data - _not_ neccesarily sufficient for restart
*/
void CRotParticle::saveSnapShotData(std::ostream& oStream)
{
  CParticle::saveSnapShotData(oStream);
  const char delim = ' ';
  oStream
    << delim
    << getQuat()
    << delim
    << getAngVel();
}



/*!
  Save check-point data sufficient for restart
*/
void CRotParticle::saveCheckPointData(std::ostream& oStream)
{
  CParticle::saveCheckPointData(oStream);
  const char delim = ' ';
  oStream
    << delim << m_inertRot
    << delim << getQuat()
    << delim << getAngVel()
    << delim << m_is_rot;
}

/*!
  load data saved with CRotParticle::saveCheckPointData

  \param iStream the input stream
*/
void CRotParticle::loadCheckPointData(std::istream &iStream)
{
  CParticle::loadCheckPointData(iStream);
  iStream >> m_inertRot;
  m_div_inertRot = 1.0/m_inertRot ;
  iStream >> m_quat;
  iStream >> m_angVel;
  iStream >> m_is_rot;
}

/*!
  get deviatoric stress
*//*
double CRotParticle::sigma_d() const
{
  double scale=1.0/(M_PI*m_rad*m_rad);
                                                                                
  double sig_d=(m_sigma-m_sigma.trace()*Matrix3::Unit()).norm();
                                                                                
  return sig_d;
}*/

ostream& operator<<(ostream& ost, const CRotParticle& CP)
{
  ost << "--CParticle " << CP.m_global_id << "  --\n";
  ost << "Radius : " << CP.m_rad << " Mass : " << CP.m_mass << "\n";
  ost << "Position : " << CP.m_pos << "\n";
  ost << "Velocity : " << CP.m_vel << "\n";
  ost << "Force    : " << CP.m_force << "\n";
  return ost;
}

