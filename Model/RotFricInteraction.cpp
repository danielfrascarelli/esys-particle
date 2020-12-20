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

#include "Model/RotFricInteraction.h"
#include "Foundation/console.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

extern double calc_angle(double , double);

CRotFrictionIGP::CRotFrictionIGP()
  : AIGParam(),
    k(0.0),
    mu_d(0.0),
    mu_s(0.0),
    k_s(0.0),
    dt(0.0),
    scaling(true),
    rigid(false),
    meanR_scaling(true)
{
}

CRotFrictionIGP::CRotFrictionIGP(
  const std::string &name,
  double kR,
  double muD,
  double muS,
  double kS,
  double dT,
  bool   scaling,
  bool   rigid,
  bool   meanR_scaling
)
  : AIGParam(name),
    k(kR),
    mu_d(muD),
    mu_s(muS),
    k_s(kS),
    dt(dT),
    scaling(scaling),
    rigid(rigid),
    meanR_scaling(meanR_scaling)
{
}

CRotFrictionIGP::CRotFrictionIGP(
  const  std::string &name,
  double youngsModulus,
  double poissonsRatio,
  double mu_d,
  double mu_s,
  double dt,
  bool   rigid,
  bool   meanR_scaling
)
  : AIGParam(name),
    mu_d(mu_d),
    mu_s(mu_s),
    dt(dt),
    scaling(true),
    rigid(rigid),
    meanR_scaling(meanR_scaling)
{
    double shearModulus = youngsModulus / (2.0*(1. + poissonsRatio));

    k = M_PI*youngsModulus;
    k_s = M_PI*shearModulus;
}

void CRotFrictionIGP::setTimeStepSize(double timeStepSize)
{
  this->dt = timeStepSize;
}

CRotFrictionInteraction::CRotFrictionInteraction():ARotPairInteraction() 
{
  m_k=0.0;
  m_mu_d=0.0;
  m_mu_s=0.0;
  m_r0=0.0;
  m_ks=0.0;
  m_dt=0.0;
  m_is_slipping=false;
  m_is_touching=false;
  m_Ffric = Vec3(0.0,0.0,0.0);
  m_E_diss = 0.0; 
  m_scaling = true;
  m_meanR_scaling = true;
  m_rigid = false;
}

CRotFrictionInteraction::CRotFrictionInteraction(
  CRotParticle* p1,
  CRotParticle* p2,
  const CRotFrictionIGP& param
)
  : ARotPairInteraction(p1,p2)
{
  double effR=1.0;
  double effL=1.0;
  double effA=1.0;


  m_r0=p1->getRad()+p2->getRad();

  m_scaling = param.scaling;
  m_meanR_scaling = param.meanR_scaling;
  m_rigid = param.rigid;
  // scale elastic param
  if (m_scaling == true) {
    if(!CParticle::getDo2dCalculations()){
      if (m_meanR_scaling == true) {
        effR=0.5*m_r0;
      }
      else {
        effR=fmin(p1->getRad(), p2->getRad());
      }
      effL=m_r0;
      effA = effR * effR;
    }
  }
  m_k =  param.k * effA / effL;
  m_ks =  param.k_s * effA / effL;
 
  m_mu_d = param.mu_d;
  m_mu_s = param.mu_s;
 
  m_dt=param.dt;
  m_cpos=p1->getPos()+((p2->getPos()-p1->getPos())*p1->getRad()/m_r0);
  m_is_slipping=false;
  m_is_touching=false;
  
  m_Ffric = Vec3(0.0,0.0,0.0);
  m_E_diss = 0.0; 
}

void CRotFrictionInteraction::setTimeStepSize(double dt)
{
  m_dt = dt;
}

CRotFrictionInteraction::~CRotFrictionInteraction()
{
}


void CRotFrictionInteraction::calcForces()
{
  if (m_rigid) {
     calcRigidBodyForces() ;
  }
  else {
     calcSimpleForces() ;
  }
} 
/*!
  Calculate elastic and frictional forces. 
  Rigid body rotations of particle-pairs are ignored
*/
void CRotFrictionInteraction::calcSimpleForces()
{
  // set current coefficient of friction depending on static/dynamic status
  double mu_current= m_is_slipping ? m_mu_d : m_mu_s ;
  // calculate distance
  const Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  const double eq_dist=m_p1->getRad()+m_p2->getRad();
  // check if there is contact
  if(dist<(eq_dist*eq_dist)){ // contact -> calculate forces
    //--- elastic force ---
    dist=sqrt(dist);
    const Vec3 force(D*(m_k*(dist-eq_dist)/dist));
    m_normal_force=force;
    const Vec3 pos(m_p2->getPos()+(m_p2->getRad()/eq_dist)*D);

    // apply elastic force 
    m_p2->applyForce(force,pos);
    m_p1->applyForce(-1.0*force,pos); 

    //--- frictional force ---
    // Calculate the relative linear-displacement of point
    // pos on each of the particles.
    const Vec3 vp1 = m_p1->getVel() + cross(m_p1->getAngVel(), pos-m_p1->getPos());
    const Vec3 vp2 = m_p2->getVel() + cross(m_p2->getAngVel(), pos-m_p2->getPos());
    Vec3 ds((vp2-vp1)*m_dt);

    // Get the tangential component by subtracting off the normal component.
    ds -= ((ds*D)/(D.norm2()))*D;

    m_Ffric += m_ks*ds;
    const double FfricNorm = m_Ffric.norm();    
    const double forceNorm = force.norm();
    if (FfricNorm > forceNorm*mu_current) { // tangential force greater than current friction -> dynamic
      m_Ffric = m_Ffric*((m_mu_d*forceNorm)/FfricNorm);
      m_force_deficit = Vec3(0.0,0.0,0.0);
      m_is_slipping = true;
      m_E_diss = fabs(m_Ffric * ds); // energy dissipated
    } else { // static friction or no frictional force
      m_is_slipping = false;
      m_E_diss = 0.0; // no energy dissipated
    }

    const Vec3 Moment1(cross(pos-m_p1->getPos(),  m_Ffric));
    const Vec3 Moment2(cross(pos-m_p2->getPos(), -m_Ffric));
    m_p1->applyMoment(Moment1);
    m_p2->applyMoment(Moment2);

    m_p1->applyForce(m_Ffric,pos);
    m_p2->applyForce(-1.0*m_Ffric,pos);
    m_cpos=pos;
    m_is_touching=true;
  } else { // no contact -> all forces are 0
    m_Ffric=Vec3(0.0,0.0,0.0);
    m_force_deficit=Vec3(0.0,0.0,0.0);
    m_normal_force=Vec3(0.0,0.0,0.0);
    m_is_slipping=false;
    m_is_touching=false;
    m_E_diss=0.0; // no energy dissipated
  }
}

/**
 * Yucang Wang's friction implementation which takes into
 * account rigid body rotation of particle-pairs.
 */
void CRotFrictionInteraction::calcRigidBodyForces()
{  //cout << "wyc in RotFric:: calcF " <<endl;
  Vec3 pos;
  Vec3 force;
  Vec3 dv,ds ;
  Vec3  d_Ffric ; 
  // calculate distance
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  // check if there is contact
  if(dist<(eq_dist*eq_dist)){ // contact -> calculate forces
    //--- elastic force ---
    dist=sqrt(dist);
    force=D*(m_k*(dist-eq_dist)/dist);
    m_normal_force=force;
    pos=m_p2->getPos()+(m_p2->getRad()/eq_dist)*D;

    // apply elastic force
    m_p1->applyForce(-1.0*force,pos);
    m_p2->applyForce(force,pos);

    //--- frictional force ---

    const Vec3 vp1_trs = m_p1->getVel();
    const Vec3 vp2_trs = m_p2->getVel();
    const Vec3 vp1_rot = cross(m_p1->getAngVel(),pos-m_p1->getPos());
    const Vec3 vp2_rot = cross(m_p2->getAngVel(),pos-m_p2->getPos());

    const Vec3 dv_trs =   vp2_trs-vp1_trs ;
    //   tangential part
    const Vec3 dv_trs_s = dv_trs - (dot(dv_trs,D)/D.norm2())*D ;
    const Vec3 ds = (vp2_rot-vp1_rot + dv_trs_s )*m_dt   ;

    if(!m_is_slipping)
    {
      if (!m_is_touching)
      {
        m_Ffric = Vec3(0.0,0.0,0.0); //first touch
      }
      else
      {

        //         due to  motion of 2 particles as a rigid body !
        //         Matrix3 mat0 = (m_p1->getQuat()).to_matrix() ;
        Vec3 rbp  = m_p2->getPos() - m_p1->getPos() ;
        Vec3 vbp  = m_p2->getVel() - m_p1->getVel() ;
        double rbp0 = rbp.norm() ;

        Vec3 omiga_s = 0.5*( m_p1->getAngVel()+m_p2->getAngVel());
        Vec3 omiga_spin  = dot(omiga_s,rbp)*rbp/(rbp0*rbp0) ;
        Vec3 omiga_m = cross(rbp,vbp)/(rbp0*rbp0) ;

        d_Ffric = m_dt*cross(omiga_spin + omiga_m , m_Ffric )  ;
        m_Ffric += d_Ffric ;

      }

      if((m_Ffric+m_ks*ds).norm()>force.norm()*m_mu_s)
      {
        m_is_slipping=true ;
        m_Ffric= m_mu_d*force.norm() *ds/ds.norm();
        m_force_deficit = Vec3(0.0,0.0,0.0);
        m_E_diss = fabs(m_Ffric * ds); // energy dissipated
      }
      else
      {
        m_Ffric += m_ks*ds;
        m_E_diss = 0.0; // no energy dissipated
      }
    }
    else
    {
      if (ds.norm() > 1.0e-8)
      {
        m_Ffric= m_mu_d*force.norm() *ds/ds.norm();
        m_force_deficit = Vec3(0.0,0.0,0.0);
        m_E_diss = fabs(m_Ffric * ds); // energy dissipated
      }
      else
      {
        m_is_slipping= false;
        m_E_diss = 0.0; // no energy dissipated
      }
    }

    const Vec3 Moment1(cross(pos-m_p1->getPos(),  m_Ffric));
    const Vec3 Moment2(cross(pos-m_p2->getPos(), -m_Ffric));

    m_p1->applyForce(m_Ffric,pos);
    m_p2->applyForce(-1.0*m_Ffric,pos);
    m_p1->applyMoment(Moment1) ;
    m_p2->applyMoment(Moment2) ;
    m_cpos=pos;
    m_is_touching=true;
  }
  else
  {
    // no contact -> all forces are 0
    m_Ffric=Vec3(0.0,0.0,0.0);
    m_force_deficit=Vec3(0.0,0.0,0.0);
    m_normal_force=Vec3(0.0,0.0,0.0);
    m_is_slipping=false;
    m_is_touching=false;
    m_E_diss = 0.0; // no energy dissipated
  }
}

bool CRotFrictionInteraction::isPersistent()
{
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  return dist<=(eq_dist*eq_dist);
}

/*!
  get the force needed to overcome friction and make the interaction slip
*/
double CRotFrictionInteraction::getAbsForceDeficit()const
{
  return m_force_deficit.norm();
}

/*!
 Calculate the normal force.
*/
void CRotFrictionInteraction::calcNormalForce()
{
  Vec3 pos;
  // calculate distance
  const Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  const double eq_dist=m_p1->getRad()+m_p2->getRad();
  // check if there is contact
  if(dist<(eq_dist*eq_dist)){ // contact -> calculate forces
    //--- elastic force ---
    dist=sqrt(dist);
    m_normal_force=D*(m_k*(dist-eq_dist)/dist);
  }
}

/*!
  get the potential energy stored in the interaction

  \warning For performance reasons the tangential part of the elastic energy is calculated directly from the current tangential force whereas it would be more accurate to calculate it incrementally during the force calculation. Be aware that therefore the potential energy is an approximation. Tests suggest an accuracy of a few percent.
*/
double CRotFrictionInteraction::getPotentialEnergy() const
{
  double e_pot_norm=0.5*m_normal_force*m_normal_force/m_k;
  double e_pot_tan=0.5*m_Ffric*m_Ffric/m_ks;

  return e_pot_norm+e_pot_tan;
}

/*!
  Get the static/dynamic status of the interaction. Returns 1 for a contact in dynamic 
  friction, 0 for static or no contact
*/
double CRotFrictionInteraction::getSlipping() const
{
  double res=m_is_slipping ? 1.0 : 0.0;
  return res;
}

/*!
  Get "sticking" contacts, i.e. return 1 if the contact is touching but not 
  slipping, 0 otherwise
*/
double CRotFrictionInteraction::getSticking() const
{
  const double res=(m_is_touching && !m_is_slipping) ? 1.0 : 0.0;
  return res;
}

/*!
  return the amount of energy dissipated during the last time step
*/
double CRotFrictionInteraction::getDissipatedEnergy() const
{
  return m_E_diss;
}

/*!
  If the particles are in contact, get total force, if not in contact return (0,0,0)
*/
Vec3 CRotFrictionInteraction::getForce() const
{
  const Vec3 f=m_is_touching ? m_Ffric-m_normal_force : Vec3(0.0,0.0,0.0);
  return f; 
}

/*!
  If the particles are in contact, get normal force, if not in contact return (0,0,0)
*/
Vec3 CRotFrictionInteraction::getNormalForce() const
{
  const Vec3 f=m_is_touching ? m_normal_force : Vec3(0.0,0.0,0.0);
  return f; 
}

Vec3 CRotFrictionInteraction::getTangentialForce() const
{
  const Vec3 f=m_is_touching ? m_Ffric : Vec3(0.0,0.0,0.0);
  return f; 
}

/*!
  return 1 if particles are in contact, 0 otherwise
*/
double CRotFrictionInteraction::Count() const
{
  double res=m_is_touching ? 1.0 : 0.0;

  return res;
}

/*!
  Return distance slipped at the contact in the current time step. 
  Returns 0 if interaction is sticking or if particles are not in contact.
*/
double CRotFrictionInteraction::getAbsSlip() const
{
        double slipdist=0.0;
        
        if(m_is_touching && m_is_slipping){
                // calculate distance
                const Vec3 D=m_p1->getPos()-m_p2->getPos();
                
                // Calculate the relative linear-displacement of point
                // pos on each of the particles.
                const Vec3 vp1 = m_p1->getVel() + cross(m_p1->getAngVel(), m_cpos-m_p1->getPos());
                const Vec3 vp2 = m_p2->getVel() + cross(m_p2->getAngVel(), m_cpos-m_p2->getPos());
                Vec3 ds((vp2-vp1)*m_dt);

                // Get the tangential component by subtracting off the normal component.
                ds -= ((ds*D)/(D.norm2()))*D;
                
                // Get slip distance from tangential velocity and time step
                slipdist=ds.norm()*m_dt;
        }
        return slipdist;
}

/*!
  Get the particle member function which returns a scalar field of a given name.

  \param name the name of the field 
*/
CRotFrictionInteraction::ScalarFieldFunction CRotFrictionInteraction::getScalarFieldFunction(const string& name)
{
  CRotFrictionInteraction::ScalarFieldFunction sf;

  if(name=="force_deficit"){
    sf=&CRotFrictionInteraction::getAbsForceDeficit;
  } else if (name=="potential_energy"){
    sf=&CRotFrictionInteraction::getPotentialEnergy;
  } else if (name=="slipping"){
    sf=&CRotFrictionInteraction::getSlipping;
  } else if (name=="sticking"){
    sf=&CRotFrictionInteraction::getSticking;
  } else if (name=="count"){
    sf=&CRotFrictionInteraction::Count;
  } else if (name=="dissipated_energy") {
    sf=&CRotFrictionInteraction::getDissipatedEnergy;
  } else if (name=="slip_distance") {
     sf=&CRotFrictionInteraction::getAbsSlip;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction scalar access function" << endl; 
  }
  
  return sf;
}

CRotFrictionInteraction::CheckedScalarFieldFunction CRotFrictionInteraction::getCheckedScalarFieldFunction(const string& name)
{
  CRotFrictionInteraction::CheckedScalarFieldFunction sf = NULL;
  cerr << "ERROR - invalid name for interaction scalar access function" << endl; 
  return sf;
}

/*!
  Get the particle member function which returns a vector field of a given name.

  \param name the name of the field 
*/
CRotFrictionInteraction::VectorFieldFunction CRotFrictionInteraction::getVectorFieldFunction(const string& name)
{
  CRotFrictionInteraction::VectorFieldFunction vf=NULL;

  if (name=="force") {
    vf = &CRotFrictionInteraction::getForce;    
  } else if (name=="normal_force") {
    vf = &CRotFrictionInteraction::getNormalForce;    
  } else if (name=="tangential_force") {
    vf = &CRotFrictionInteraction::getTangentialForce;    
  } else if (name=="frictional_force") {
    vf = &CRotFrictionInteraction::getTangentialForce;    
  } else {
    cerr << "ERROR - invalid name for interaction vector access function" << endl;
  }

  return vf;
}

/*!
  Pack a CFrictionInteraction into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<CRotFrictionInteraction>(const CRotFrictionInteraction& I)
{
  append(I.m_k);
  append(I.m_r0);
  append(I.m_mu_d);
  append(I.m_mu_s);
  append(I.m_ks);
  append(I.m_dt);
  append(static_cast<int>(I.m_scaling));
  append(static_cast<int>(I.m_meanR_scaling));
  append(I.m_id[0]);
  append(I.m_id[1]);
  append(static_cast<int>(I.m_is_slipping));
  append(I.m_Ffric.X());
  append(I.m_Ffric.Y());
  append(I.m_Ffric.Z());
}

/*!
  Unpack a CFrictionInteraction from a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<CRotFrictionInteraction>(CRotFrictionInteraction& I)
{
  I.m_k=pop_double();
  I.m_r0=pop_double();
  I.m_mu_d=pop_double();
  I.m_mu_s=pop_double();
  I.m_ks=pop_double();
  I.m_dt=pop_double();
  I.m_scaling=static_cast<bool>(pop_int());
  I.m_meanR_scaling=static_cast<bool>(pop_int());
  I.m_id.clear();
  I.m_id.push_back(pop_int());
  I.m_id.push_back(pop_int());
  I.m_is_slipping = static_cast<bool>(pop_int());
  I.m_Ffric.X() = pop_double();
  I.m_Ffric.Y() = pop_double();
  I.m_Ffric.Z() = pop_double();
}

/*!
  Save restart data to an open ostream 

  \param oStream the output stream
*/
void CRotFrictionInteraction::saveRestartData(std::ostream &oStream)
{
  oStream << m_k << " ";
  oStream << m_r0 << " ";
  oStream << m_mu_d << " ";
  oStream << m_mu_s << " ";
  oStream << m_ks << " ";
  oStream << m_dt << " ";
  oStream << m_scaling << " ";
  oStream << m_meanR_scaling << " ";
  oStream << m_id[0] << " ";
  oStream << m_id[1] << " ";
  oStream << m_is_slipping << " ";
  oStream << m_is_touching << " ";
  oStream << m_Ffric.X() << " ";
  oStream << m_Ffric.Y() << " ";
  oStream << m_Ffric.Z();
}


/*!
  Load restart data from an open istream 

  \param iStream the input stream
*/
void CRotFrictionInteraction::loadRestartData(std::istream &iStream)
{
  iStream >> m_k ;
  iStream >> m_r0 ;
  iStream >> m_mu_d ;
  iStream >> m_mu_s ;
  iStream >> m_ks ;
  iStream >> m_dt ;
  iStream >> m_scaling ;
  iStream >> m_meanR_scaling ;
  iStream >> m_id[0] ;
  iStream >> m_id[1] ;
  iStream >> m_is_slipping ;
  iStream >> m_is_touching ;
  iStream >> m_Ffric.X() ;
  iStream >> m_Ffric.Y() ;
  iStream >> m_Ffric.Z();
}

ostream& operator<<(ostream& ost,const CRotFrictionInteraction& BI)
{
  ost << "[" << BI.m_p1->getID() << " - " << BI.m_p2->getID() << "]";
  return ost;
}
