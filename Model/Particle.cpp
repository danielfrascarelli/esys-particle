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


#include "Parallel/mpibuf.h"
#include "tml/message/packed_message_interface.h"
#include "Parallel/mpisgbuf.h"

#include "Model/Particle.h"
#include "Geometry/SimpleParticleData.h"
#include "Foundation/console.h"

const Vec3 zeroVec3(0.0,0.0,0.0);

bool CParticle::s_do2Calculations = false;

CParticle::CParticle()
  : CBasicParticle(zeroVec3,0.0),
    m_sigma(), 
    m_vel(),
    m_force(),
    m_oldpos(),
    m_initpos(),
    m_circular_shift(), //!< shift vector if particle is circular image
    m_mass(0.0),
    m_div_mass(0.0),
    m_div_vol(0.0),
    flag(false),
  m_is_dynamic(true)
{
}

CParticle::CParticle(const esys::lsm::SimpleParticleData &data)
  : CBasicParticle(data),
    m_sigma(),
    m_vel(),
    m_force(),
    m_oldpos(data.getPosition()),
    m_initpos(data.getPosition()),
    m_circular_shift(), //!< shift vector if particle is circular image
    m_mass(data.getMass()),
    m_div_mass((data.getMass() != 0) ? 1.0/data.getMass() : 0.0),
    m_div_vol(1.0/(4.0/3.0*M_PI*data.getRadius()*data.getRadius()*data.getRadius())),
    flag(false),
    m_is_dynamic(true) //! SimpleParticleData has no is_dynamic info -> default to true 
{}

/*!
  Construct particle. Old and initial position are assumed to be identical to current position. 

  \param rad radius
  \param mass mass
  \param pos current position
  \param vel current velocity
  \param force currently applied force
  \param id particle id
*/
CParticle::CParticle(double rad,double mass,const Vec3& pos,const Vec3& vel,const Vec3& force,int id,bool is_dyn) : CBasicParticle(pos,rad)
{
  m_oldpos=pos;
  m_initpos=pos;
  m_vel=vel;
  m_force=force;
  m_mass=mass;
  if(m_mass!=0.0){
    m_div_mass=1.0/m_mass;
  } else {
    m_div_mass=0.0;
  }
  m_div_vol=1.0/(4.0/3.0*M_PI*rad*rad*rad);
  m_global_id=id;
  m_circular_shift = zeroVec3;
  flag=false;
  m_is_dynamic=is_dyn;
}

CParticle::CParticle(double rad,double mass,const Vec3& pos,const Vec3& oldpos,const Vec3& initpos,const Vec3& vel,const Vec3& force,int id,bool is_dyn) : CBasicParticle(pos,rad)
{
  m_oldpos=oldpos;
  m_initpos=initpos;
  m_vel=vel;
  m_force=force;
  m_mass=mass;
  if(m_mass!=0.0){
    m_div_mass=1.0/m_mass;
  } else {
    m_div_mass=0.0;
  }
  m_div_vol=1.0/(4.0/3.0*M_PI*rad*rad*rad);
  m_global_id=id;
  m_circular_shift = zeroVec3;
  flag=false;
  m_is_dynamic=is_dyn;
}

/*!
  set the density of the particle
*/
void CParticle::setDensity(double rho)
{
  const double pi=3.141592654;
  
  if(getDo2dCalculations()){
    m_mass=rho*pi*m_rad*m_rad;
  } else {
    m_mass=(4.0/3.0)*rho*pi*m_rad*m_rad*m_rad;
  }
  if(m_mass!=0.0){
    m_div_mass=1.0/m_mass;
  } else {
    m_div_mass=0.0;
  }
}

/*!
  Get the particle member function which returns a scalar field of a given name.

  \param name the name of the field 
*/
CParticle::ScalarFieldFunction CParticle::getScalarFieldFunction(const string& name)
{
  CParticle::ScalarFieldFunction sf;

  if(name=="id"){
    sf=&CParticle::getIDField;
  } else if(name=="tag"){
    sf=&CParticle::getTagField;    
  } else if(name=="sigma_xx_2d"){
    sf=&CParticle::sigma_xx_2D;
  } else if(name=="sigma_xy_2d"){
    sf=&CParticle::sigma_xy_2D;
  } else if(name=="sigma_yy_2d"){
    sf=&CParticle::sigma_yy_2D;
  } else if(name=="sigma_d"){
    sf=&CParticle::sigma_d;
  } else if(name=="sigma_vonMises"){
    sf=&CParticle::sigma_vonMises;
  } else if(name=="e_kin"){
    sf=&CParticle::getKineticEnergy;
  } else if(name=="radius"){
    sf=&CParticle::getRad;
  } else if(name=="v_abs"){
    sf=&CParticle::getAbsVel;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for particle scalar  access function" << endl; 
  }

  return sf;
}

/*!
  Get the particle member function which returns a vector field of a given name.

  \param name the name of the field 
*/
CParticle::VectorFieldFunction CParticle::getVectorFieldFunction(const string& name)
{
  CParticle::VectorFieldFunction sf;

  if(name=="displacement"){
    sf=&CParticle::getTotalDisplacement;
  } else if(name=="velocity"){
    sf=&CParticle::getVel;
  } else if(name=="position"){
    sf=&CParticle::getPos;
  } else if(name=="force"){
    sf=&CParticle::getForce;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for particle vector access function" << endl; 
  }

  return sf;
}

/*!
  Pack a CParticle into a TML packed message

  \param p the particle
  \todo BasicParticle data should be handled by pack<Basicparticle>
*/
template<>
void TML_PackedMessageInterface::pack<CParticle>(const CParticle& p)
{

  append(p.m_tag);
  append(p.m_global_id);
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
  append(p.m_is_dynamic);
  append(p.m_div_vol);
}

/*!
  Unpack a CParticle from a TML packed message

  \param p the particle
*/
template<>
void TML_PackedMessageInterface::unpack<CParticle>(CParticle& p)
{
  double db[6*3+2] ;

  p.m_tag=pop_int();
  p.m_global_id=pop_int();
  pop_doubles(db,20);
  p.m_pos=Vec3(db[0],db[1],db[2]) ;
  p.m_oldpos=Vec3(db[3],db[4],db[5]) ;
  p.m_initpos=Vec3(db[6],db[7],db[8]) ;
  p.m_circular_shift=Vec3(db[9],db[10],db[11]) ;
  p.m_vel=Vec3(db[12],db[13],db[14]) ;
  p.m_force=Vec3(db[15],db[16],db[17]) ;
  p.m_rad = db[18];
  p.m_mass = db[19];
  if(p.m_mass!=0.0){
    p.m_div_mass=1.0/p.m_mass;
  } else {
    p.m_div_mass=0.0;
  }
  p.m_is_dynamic=pop_bool();
  p.m_div_vol=pop_double();
}

/*!
  Do the time integration for the particle by 1st order method.

  \param dt the time step
*/
void CParticle::integrate(double dt)
{
  if(m_is_dynamic){
    if (getDo2dCalculations()) {
      m_force = Vec3(m_force.X(), m_force.Y(), 0);
    }
#if 1
    // 1st order scheme,  13 Flops
    m_vel+=m_force*m_div_mass*dt; // 7 Flops 
    m_pos+=m_vel*dt;              // 6 Flops
#else
    // 2nd order scheme.
    const Vec3  Acc = m_force*m_div_mass ;
    m_pos +=  dt*m_vel + 0.5*dt*dt*Acc ;
    m_vel +=  dt*Acc;
#endif
  }
}

/*!
  zero forces on particle
*/
void CParticle::zeroForce()
{
  m_force=Vec3(0.0,0.0,0.0);
  m_sigma=Matrix3();
}

/*!
  write particle data as a line in openDX format into a stream(file)
  
  \param ost the output stream
  \param slid from which sublattice the particle is coming
*/
void CParticle::writeAsDXLine(ostream& ost,int slid)
{
  ost << m_pos.X() << " " << m_pos.Y() << " " << m_pos.Z() << " " ; //position
  ost << slid << " " ;
  ost << m_rad << " " << m_mass << " " ; // radius,mass
  Vec3 disp=m_pos-m_initpos; 
  ost << disp.X() << " "  << disp.Y() << " "  << disp.Z() << " " ; //displacement
  ost << m_vel.X() << " "  << m_vel.Y() << " "  << m_vel.Z() << " " ; // velocity
  ost << endl;
}

/*!
  get values to be exchanged, i.e. pos and vel
*/
CParticle::exchangeType CParticle::getExchangeValues()
{
  return 
    exchangeType(
      m_pos     - m_circular_shift,
      m_initpos - m_circular_shift,
      m_oldpos  - m_circular_shift,
      m_vel,
      m_is_dynamic
    );
}

/*!
  set pos and vel from exchangeType

  \param E the exchanged values, i.e. a pair of pos and vel 
*/
void CParticle::setExchangeValues(const exchangeType& E)
{
  m_pos     = E.m_pos     + m_circular_shift;
  m_initpos = E.m_initPos + m_circular_shift;
  m_oldpos  = E.m_oldPos  + m_circular_shift;
  m_vel     = E.m_vel;
  m_is_dynamic = E.m_is_dynamic;
}

/*!
  set circular shift vector

  \param cv the circular shift vector
*/
void CParticle::setCircular(const Vec3& cv)
{
  m_circular_shift += cv;
  m_pos            += cv;
  m_initpos        += cv;
  m_oldpos         += cv;
}

/*!
  Pack an exchangeType into a TML packed message

  \param p the exchangeType
*/
template<>
void TML_PackedMessageInterface::pack<CParticle::exchangeType>(const CParticle::exchangeType& p)
{
  append(p.m_pos.X());
  append(p.m_pos.Y());
  append(p.m_pos.Z());
  append(p.m_initPos.X());
  append(p.m_initPos.Y());
  append(p.m_initPos.Z());
  append(p.m_oldPos.X());
  append(p.m_oldPos.Y());
  append(p.m_oldPos.Z());
  append(p.m_vel.X());
  append(p.m_vel.Y());
  append(p.m_vel.Z());
  append(p.m_is_dynamic);
}

/*!
  Unpack an exchangeType from a TML packed message
*/
template<>
void TML_PackedMessageInterface::unpack<CParticle::exchangeType>(CParticle::exchangeType& p)
{
  const int numDoubles = 12;
  double db[numDoubles];

  pop_doubles(db, numDoubles);
  p.m_pos     = Vec3(db[0],  db[1],  db[2]);
  p.m_initPos = Vec3(db[3],  db[4],  db[5]);
  p.m_oldPos  = Vec3(db[6],  db[7],  db[8]);
  p.m_vel     = Vec3(db[9], db[10], db[11]);
  p.m_is_dynamic=pop_bool();
}


/*!
  Apply a force to a particle at a given position. Also updates particle stress tensor. 

  \param force the force
  \param pos the position at which the force acts on the particle
*/
void CParticle::applyForce(const Vec3& force,const Vec3& pos)
{
  m_force+=force;
  Vec3 rel_pos=pos-m_pos;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      m_sigma(i,j)-=force[i]*rel_pos[j];
    }
  }
}

/*!
  get full stress tensor
*/
Matrix3 CParticle::sigma() const
{
    return m_sigma*m_div_vol;
};

/*!
  get deviatoric stress
*/
double CParticle::sigma_d() const
{
  //double scale=1.0/(M_PI*m_rad*m_rad);
  
  double sig_d=(m_sigma-(m_sigma.trace()/3.)*Matrix3::Unit()).norm();

  return sig_d;
}

/*!
  get von Mises stress
*/
double CParticle::sigma_vonMises() const
{
  return sqrt(1.5)*sigma_d();
}

ostream& operator<<(ostream& ost, const CParticle& CP)
{
  ost << "--CParticle " << CP.m_global_id << "  --\n";
  ost << "Radius : " << CP.m_rad << " Mass : " << CP.m_mass << "\n"; 
  ost << "Position : " << CP.m_pos << "\n";
  ost << "Velocity : " << CP.m_vel << "\n";
  ost << "Force    : " << CP.m_force << "\n";
  return ost;
}

/*!
  save snapshot data _not_ neccesarily sufficient for restart
*/
void CParticle::saveSnapShotData( std::ostream& oStream){
  const char delim = ' ';
  oStream 
    << getPos()                        << delim // 0,1,2
    << getRad()                        << delim // 3
    << getID()                         << delim // 4
    << getTag()                        << delim // 5
    << getMass()                       << delim // 6
    << getInitPos() - m_circular_shift << delim // 7,8,9
    << getOldPos()                     << delim // 10,11,12
    << getVel()                        << delim // 13,14,15
    << getForce()                      << delim // 16,17,18
    << m_circular_shift;                        // 19,20,21
}

/*!
  save checkpoint data, sufficient for restart
*/
void CParticle::saveCheckPointData( std::ostream& oStream){
  const char delim = ' ';
  oStream 
    << getPos()                        << delim // 0,1,2
    << getRad()                        << delim // 3
    << getID()                         << delim // 4
    << getTag()                        << delim // 5
    << getMass()                       << delim // 6
    << getInitPos() - m_circular_shift << delim // 7,8,9
    << getOldPos()                     << delim // 10,11,12
    << getVel()                        << delim // 13,14,15
    << getForce()                      << delim // 16,17,18
    << m_circular_shift << delim                // 19,20,21
    << m_is_dynamic;
}

/*!
  load checkpoint data saved with CParticle::saveCheckPointData
*/
void CParticle::loadCheckPointData(std::istream &iStream)
{
  Vec3 pos;
  Vec3 oldpos;
  Vec3 initpos;
  Vec3 vel;
  Vec3 force;
  Vec3 circular_shift;
  double rad;
  double mass;
  int id;
  int tag;
  bool is_dyn;

  iStream
    >> pos
    >> rad
    >> id
    >> tag
    >> mass
    >> initpos
    >> oldpos
    >> vel
    >> force
    >> circular_shift
    >> is_dyn;

  CParticle 
    particle
    (
      rad,
      mass,
      pos,
      oldpos,
      initpos + circular_shift,
      vel,
      force,
      id, is_dyn
    );
  particle.setTag(tag);
  
  *this = particle;
}
