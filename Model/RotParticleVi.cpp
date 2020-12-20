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
#include "Model/RotParticleVi.h"
#include "Geometry/SimpleParticleData.h"
#include "Foundation/console.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

CRotParticleVi::CRotParticleVi() : CParticle()
{
  m_global_id=-1;
  flag=false;
}

CRotParticleVi::CRotParticleVi(const esys::lsm::SimpleParticleData &particleData) : CParticle(particleData)
{
  m_quat = Quaternion(1.0,Vec3(0.0,0.0,0.0));
  m_initquat = m_quat;
  m_inertRot  = 0.5*particleData.getMass()*particleData.getRadius()*particleData.getRadius();
  m_div_inertRot = 1.0/m_inertRot;
  m_angVel = Vec3(0.0,0.0,0.0);
  m_angVel_t = Vec3(0.0,0.0,0.0) ;
  m_moment     = Vec3(0.0,0.0,0.0);
}

CRotParticleVi::CRotParticleVi(
  double rad,
  double mass,
  const Vec3& pos,
  const Vec3& vel,
  const Vec3& force,
  int id,
  bool is_dyn
) : CParticle(rad, mass, pos, vel, force, id, is_dyn)
{
  m_quat = Quaternion(1.0,Vec3(0.0,0.0,0.0));
  m_initquat = m_quat;
  m_inertRot  = 0.5*mass*rad*rad;
  m_div_inertRot = 1.0/m_inertRot;
  m_angVel = Vec3(0.0,0.0,0.0);
  m_angVel_t = Vec3(0.0,0.0,0.0) ;
  m_moment     = Vec3(0.0,0.0,0.0);
  m_is_dynamic=is_dyn;
}

CRotParticleVi::CRotParticleVi(const CParticle &particle) : CParticle(particle)
{
  m_quat = Quaternion(1.0,Vec3(0.0,0.0,0.0));
  m_initquat = m_quat;
  m_inertRot  = 0.5*particle.getMass()*particle.getRad()*particle.getRad();
  m_div_inertRot = 1.0/m_inertRot;
  m_angVel = Vec3(0.0,0.0,0.0);
  m_moment     = Vec3(0.0,0.0,0.0);
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
CRotParticleVi::CRotParticleVi(
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
  const Vec3& angvel_t
) :CParticle(rad, mass,pos,vel,force,id,true)
{
  m_circular_shift=Vec3(0.0,0.0,0.0);
  flag=false;

  m_quat         = quat;
  m_initquat     = quat;
  m_inertRot     = inertRot;
  m_moment       = moment;
  m_angVel       = angvel;
  m_angVel_t       = angvel_t; 
  m_div_inertRot = 1.0/m_inertRot;
}

CRotParticleVi::CRotParticleVi(
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
  const Vec3& angvel_t
) :CParticle(rad, mass,pos,oldpos, initpos,vel,force,id,true)
{
  m_circular_shift=Vec3(0.0,0.0,0.0);
  flag=false;

  m_quat         = quat;
  m_initquat     = initquat;
  m_inertRot     = inertRot ;
  m_moment       = moment;
  m_angVel       = angvel;
  m_angVel_t       = angvel_t;
  m_div_inertRot = 1.0/m_inertRot ;
}

/*!
  Pack a CParticle into a TML packed message

  \param p the particle
  \todo BasicParticle data should be handled by pack<Basicparticle>
*/
template<>
void TML_PackedMessageInterface::pack<CRotParticleVi>(const CRotParticleVi& p)
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
  

  append(p.m_global_id); // original one here
}


/*!
  Unpack a CParticle from a TML packed message

  \param p the particle
*/

template<>
void TML_PackedMessageInterface::unpack<CRotParticleVi>(CRotParticleVi& p)
{
  const int numElems = 38 ;
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


  p.m_global_id=pop_int();
}

/*!
  Do the time integration for the particle.

  \param dt the time step
*/


/*
void CRotParticleVi::integrate(double dt)
{
  if (CParticle::getDo2dCalculations()) {
    m_force  = Vec3(m_force.X(), m_force.Y(), 0);
    m_moment = Vec3(0, 0, m_moment.Z());
  }

  // integrate rotational part
  m_angVel += (dt*m_div_inertRot) * m_moment;
 // m_quat   += (dt/2.0)*(Quaternion(0, m_angVel)*m_quat);
  // wrong ???

   m_quat   += (dt/2.0)*(m_quat*Quaternion(0, m_angVel));

  // Integrate linear part
  CParticle::integrate(dt);
}
*/



void CRotParticleVi::integrate(double dt)
{
// need to test this 


//zjr 
//    cout << " wyc i=  " <<getID() << endl;
//    cout << " rad= " << getRad() << endl;
//    cout << " mass=  " <<1.0/m_div_mass << "   " <<1.0/m_div_inertRot << endl;
//zjr   
//  cout <<"wyc in integ  force  " << m_force << "  moment   " << m_moment << endl;

//zjr 
//    cout <<"quat=   " << m_quat << endl; 

   // m_lastquat = m_quat ;  no use  ?
   // m_lastpos  = m_pos ;


// no rotatoin !!

//    m_div_inertRot = 0.0 ; 


    if(CParticle::getDo2dCalculations()) {
   
       m_force =  Vec3( m_force.X(), m_force.Y(), 0.0);
       m_moment = Vec3(0.0,0.0,m_moment.Z()) ;
     } 


    Vec3 ang_vel_t =  m_angVel+ 0.5 * dt * m_moment*m_div_inertRot ;
    m_angVel_t = ang_vel_t ;

//    Quaternion dq = 0.5*m_quat*Quaternion(0, ang_vel_t);
//     ang_vel_t is in space-fixed frame, if it is in body-fixed frame, the above is right !

    Quaternion dq = 0.5*Quaternion(0, ang_vel_t)*m_quat ;
    
    Vec3 ang_vel_t2 =  m_angVel +  dt * m_moment*m_div_inertRot ;

    m_angVel = ang_vel_t2 ;

    Quaternion q_t2 = m_quat + 0.5*dt*dq ;
   
    Quaternion dq_t2 = 0.5*Quaternion(0, ang_vel_t2)*q_t2  ;

     m_quat +=  dt*dq_t2 ;
 

       m_vel+=m_force*m_div_mass*dt;
       m_pos+=m_vel*dt;

//      cout << " 1/mass = " << m_div_mass << "  1/iner=  " <<m_div_inertRot<< endl;

//     Vec3  Acc = m_force*m_div_mass ;
//     m_pos +=  dt*m_vel + 0.5*dt*dt*Acc ;
//     m_vel +=  dt*Acc ;
}


/*!
  zero forces on particle
*/
void CRotParticleVi::zeroForce()
{
  m_force = Vec3(0.0,0.0,0.0);
  m_moment = Vec3(0.0,0.0,0.0);
  m_sigma = Matrix3();
}

void CRotParticleVi::rescale()
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

void CRotParticleVi::applyMoment(const Vec3& moment)
{
  m_moment += moment;
}

map<string,AField*> CRotParticleVi::generateFields(ParallelParticleArray<CRotParticleVi>* particles)
{
  map<string,AField*> res;
                                                                              
  return res;
}

CRotParticleVi::exchangeType CRotParticleVi::getExchangeValues()
{
  return 
    exchangeType(
      m_pos     - m_circular_shift,
      m_initpos - m_circular_shift,
      m_vel,
      m_angVel,
      m_angVel_t,
      m_quat
    );
}

/*!
  Set pos, vel and angular vel from exchangeType
                                                                                
  \param E the exchanged values
*/
void CRotParticleVi::setExchangeValues(const exchangeType& E)
{
  m_pos        = E.m_pos     + m_circular_shift;
  m_initpos    = E.m_initPos + m_circular_shift;
  m_vel        = E.m_vel;
  m_angVel = E.m_angVel;
  m_angVel_t = E.m_angVel_t;
  m_quat       = E.m_quat;
}

template<>
void TML_PackedMessageInterface::pack<CRotParticleVi::exchangeType>(const CRotParticleVi::exchangeType& p)
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
}

/*!
  Unpack an exchangeType from a TML packed message
*/
template<>
void TML_PackedMessageInterface::unpack<CRotParticleVi::exchangeType>(CRotParticleVi::exchangeType& p)
{
  const int numDoubles = 19;
  double db[numDoubles];

  pop_doubles(db, numDoubles);
  p.m_pos        = Vec3(db[0],db[1],db[2]);
  p.m_initPos    = Vec3(db[3],db[4],db[5]);
  p.m_vel        = Vec3(db[6],db[7],db[8]);
  p.m_angVel     = Vec3(db[9],db[10],db[11]);
  p.m_angVel_t     = Vec3(db[12],db[13],db[14]);
  p.m_quat       = Quaternion(db[15],Vec3(db[16],db[17],db[18]));
}

CRotParticleVi::ScalarFieldFunction CRotParticleVi::getScalarFieldFunction(const string& name)
{
  CRotParticleVi::ScalarFieldFunction sf;
                                                                                
  if(name=="id"){
    sf=&CRotParticleVi::getIDField;
  } else if(name=="sigma_xx_2d"){
    sf=&CRotParticleVi::sigma_xx_2D;
  } else if(name=="sigma_xy_2d"){
    sf=&CRotParticleVi::sigma_xy_2D;
  } else if(name=="sigma_yy_2d"){
    sf=&CRotParticleVi::sigma_yy_2D;
  } else if(name=="sigma_d"){
    sf=&CRotParticleVi::sigma_d;
  } else if(name=="sigma_vonMises"){
    sf=&CRotParticleVi::sigma_vonMises;
  } else if(name=="e_kin"){
    sf=&CRotParticleVi::getKineticEnergy;
  } else if(name=="radius"){
    sf=&CRotParticleVi::getRad;

  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for particle scalar  access function" << endl;                                                                                
  }

  return sf;
}


CRotParticleVi::VectorFieldFunction CRotParticleVi::getVectorFieldFunction(const string& name)
{
  CRotParticleVi::VectorFieldFunction sf;
                                                                                
  if(name=="displacement"){
    sf=&CRotParticleVi::getTotalDisplacement;
  } else if(name=="velocity"){
    sf=&CRotParticleVi::getVel;
  } else if(name=="ang_velocity"){
    sf=&CRotParticleVi::getAngVelNR;
  } else if(name=="position"){
    sf=&CRotParticleVi::getPos;
  } else if(name=="force"){
    sf=&CRotParticleVi::getForce;
  } else if(name=="anglevector"){
    sf=&CRotParticleVi::getAngVector;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for particle vector access function" << endl;
  }
                                                                                
  return sf;
}



// used in 2D to output rotated angle

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
void CRotParticleVi::setCircular(const Vec3& cv)
{
  CParticle::setCircular(cv);
}

/*!
  Save snapshot data - _not_ neccesarily sufficient for restart
*/
void CRotParticleVi::saveSnapShotData(std::ostream& oStream)
{
  CParticle::saveSnapShotData(oStream);
  const char delim = ' ';
  oStream
    << delim
    << getQuat()
    << delim
    << getAngVel();
}

/**
 * Save check-point data.
 */
void CRotParticleVi::saveCheckPointData(std::ostream& oStream)
{
  CParticle::saveCheckPointData(oStream);
  const char delim = ' ';
  oStream
    << delim
    << getQuat()
    << delim
    << getAngVel();
}

void CRotParticleVi::loadCheckPointData(std::istream &iStream)
{
  CParticle::loadCheckPointData(iStream);
  iStream >> m_quat;
  iStream >> m_angVel_t;
}

/*!
  get deviatoric stress
*//*
double CRotParticleVi::sigma_d() const
{
  double scale=1.0/(M_PI*m_rad*m_rad);
                                                                                
  double sig_d=(m_sigma-m_sigma.trace()*Matrix3::Unit()).norm();
                                                                                
  return sig_d;
}*/

ostream& operator<<(ostream& ost, const CRotParticleVi& CP)
{
  ost << "--CParticle " << CP.m_global_id << "  --\n";
  ost << "Radius : " << CP.m_rad << " Mass : " << CP.m_mass << "\n";
  ost << "Position : " << CP.m_pos << "\n";
  ost << "Velocity : " << CP.m_vel << "\n";
  ost << "Force    : " << CP.m_force << "\n";
  return ost;
}

