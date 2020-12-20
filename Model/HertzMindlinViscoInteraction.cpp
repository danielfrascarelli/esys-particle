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


#include "Model/HertzMindlinViscoInteraction.h"
#include "Foundation/console.h"
#include "tml/message/packed_message_interface.h"
#include <math.h>

CHertzMindlinViscoIGP::CHertzMindlinViscoIGP()
  : AIGParam(), m_E(0.0), m_nu(0.0), mu(0.0), m_COR(0.0), dt(0.0)
{
}

CHertzMindlinViscoIGP::CHertzMindlinViscoIGP(
  const std::string &name,
  double E,
  double nu,
  double fricCoef,
  double restiCoef,
  double dT
)
  : AIGParam(name),
    m_E(E),
    m_nu(nu),
    mu(fricCoef),
    m_COR(restiCoef),
    dt(dT)
{
}

CHertzMindlinViscoInteraction
::CHertzMindlinViscoInteraction()
{
  m_E=0.0;
  m_nu=0.0;
  m_dn=0.0;
  m_mu=0.0;
  m_COR=0.0;
  m_r0=0.0;
  m_dt=0.0;
  m_is_slipping=false;
  m_is_touching=false;
  m_E_diss=0.0;
  m_E_visc_normal=0.0;
  m_E_visc_shear=0.0;
}

void CHertzMindlinViscoIGP::setTimeStepSize(double timeStepSize)
{
  this->dt = timeStepSize;
}

/*!
  constructor for CHertzMindlinViscoInteraction without friction
  parameters, only calls the constructor of APairInteraction
  with the 2 particle pointers
*/
CHertzMindlinViscoInteraction
::CHertzMindlinViscoInteraction(
  CParticle* p1,
  CParticle* p2
)
  : APairInteraction(p1,p2)
{
  m_is_slipping=false;
  m_is_touching=false;
  m_E_diss=0.0;
  m_E_visc_normal=0.0;
  m_E_visc_shear=0.0;
}


CHertzMindlinViscoInteraction
::CHertzMindlinViscoInteraction(
  CParticle* p1,
  CParticle* p2,
  const CHertzMindlinViscoIGP& param
)
  : APairInteraction(p1,p2)
{
  m_E=param.m_E;
  m_nu=param.m_nu;
  m_dn=0.0;
  m_mu=param.mu;
  m_COR=param.m_COR;
  m_r0=p1->getRad()+p2->getRad();
  m_dt=param.dt;
  m_cpos=p1->getPos()+((p2->getPos()-p1->getPos())*p1->getRad()/m_r0);
  m_is_slipping=false;
  m_is_touching=false;
  m_E_diss=0.0;
  m_E_visc_normal=0.0;
  m_E_visc_shear=0.0;
}

CHertzMindlinViscoInteraction
::~CHertzMindlinViscoInteraction()
{
}

void CHertzMindlinViscoInteraction::setTimeStepSize(double dt)
{
  m_dt = dt;
}

/*!
  Calculate viscoelastic and frictional forces.
*/
void CHertzMindlinViscoInteraction::calcForces()
{
  m_E_visc_normal = 0.0;
  m_E_visc_shear = 0.0;
  // calculate distance
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  // check if there is contact
  if(dist<(eq_dist*eq_dist))
  { // contact -> calculate forces
    //--- viscoelastic force ---
    double R_ij=1.0/(1.0/m_p1->getRad()+1.0/m_p2->getRad());
    dist=sqrt(dist);
    m_dn=eq_dist-dist;
    Vec3 dir=D.unit();
    double norm_m_normal_force=4.0/3.0*m_E*sqrt(R_ij)*pow(m_dn,1.5); //see (Renzo,2005) eq.5
    m_normal_force = dir*norm_m_normal_force;
    double M_ij=1.0/(1.0/m_p1->getMass()+1.0/m_p2->getMass()); //equivalent mass
    Vec3 v_rel=m_p2->getVel()-m_p1->getVel(); //relative velocity
    Vec3 v_nor=(v_rel*dir)*dir; //normal relative velocity
    double lnCOR= log(m_COR);
    //see (Antypov,2011) eq.6:
    Vec3 normal_damping = 2.0*sqrt(5.0/3.0)*lnCOR/sqrt(pow(lnCOR,2)+pow(3.1415926,2))*sqrt(M_ij*m_E)*pow(R_ij,0.25)*pow(m_dn,0.25)*v_nor;

    Vec3 pos=m_p2->getPos()+(m_p2->getRad()/eq_dist)*D;
    // apply elastic force and normal damping force
    m_p1->applyForce(m_normal_force-normal_damping,pos);
    m_p2->applyForce(-1.0*(m_normal_force-normal_damping),pos);
    m_E_visc_normal = (normal_damping*v_nor)*m_dt;

    //--- frictional force ---
    // particle movement since last timestep
    const Vec3 d1=m_p1->getVel()*m_dt;
    const Vec3 d2=m_p2->getVel()*m_dt;
    Vec3 ds = d2 - d1;
    // Compute tangential part by subtracting off normal component.
    //ds -= ((ds*D)/(D.norm2()))*D;
    ds -= (ds*dir)*dir;
    m_Ffric = 2.0/3.0*8.0*m_E/2.0/(1.0+m_nu)*sqrt(R_ij)*sqrt(m_dn)*ds; //combination of (Renzo,2005) eq.6 & 15, and (Tsuji, 1992) eq.20
    Vec3 v_tan = v_rel-v_nor; //tangential relative velocity
    // Compute tangential damping force:
    Vec3 tange_damping = sqrt(80.0/3.0)*lnCOR/sqrt(pow(lnCOR,2)+pow(3.1415926,2))*sqrt(M_ij*m_E/2.0/(1.0+m_nu))*pow(R_ij,0.25)*pow(m_dn,0.25)*v_tan;
    m_Ffric -= tange_damping;

    m_E_visc_shear = (tange_damping*v_tan)*m_dt;

    const double FfricNorm = m_Ffric.norm();
    const double forceNorm = m_normal_force.norm();
    // decide static/dynamic
    if (FfricNorm > forceNorm*m_mu)
    { // tangential force greater than static friction -> dynamic
      m_Ffric=m_Ffric*((m_mu*forceNorm)/FfricNorm);
      m_force_deficit=Vec3(0.0,0.0,0.0);
      m_is_slipping=true;
      m_E_diss=m_mu*fabs(m_normal_force*(d2-d1)); // energy dissipated
    }
    else if (FfricNorm > 0.0)
    { // static friction
      m_is_slipping=false;
      m_E_diss=0.0; // no energy dissipated
    }
    else
    { // no frictional force -> force deficit=mu*F_n
      m_is_slipping=false;
      m_E_diss=0.0; // no energy dissipated
    }
    m_p1->applyForce(m_Ffric, pos);
    m_p2->applyForce(-1.0*m_Ffric, pos);
    m_cpos=pos;
    m_is_touching=true;
  }
  else
  { // no contact -> all forces are 0
    m_Ffric=Vec3(0.0,0.0,0.0);
    m_normal_force=Vec3(0.0,0.0,0.0);
    m_is_slipping=false;
    m_is_touching=false;
    m_E_diss=0.0; // no energy dissipated
  }
}

bool CHertzMindlinViscoInteraction::isPersistent()
{
  const Vec3 D=m_p1->getPos()-m_p2->getPos();
  const double dist=D*D;
  const double eq_dist=m_p1->getRad()+m_p2->getRad();
  return (dist <= (eq_dist*eq_dist));
}

/*!
  get the force needed to overcome friction and make the interaction slip
*/
double CHertzMindlinViscoInteraction::getAbsForceDeficit()const
{
  return m_force_deficit.norm();
}

/*!
  get current frictional/stopping force
*/
pair<bool,double>
CHertzMindlinViscoInteraction::getAbsFrictionalForce() const
{
  pair<bool,double> res;

  // calculate distance
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  // check if there is contact
  if(dist<(eq_dist*eq_dist))
  {
    res.first=true;
    res.second=m_Ffric.norm();
  }
  else
  {
    res.first=false;
  }

  return res;
}

/*!
  get current frictional/stopping stress (f_fric/r^2)
*/
pair<bool,double>
CHertzMindlinViscoInteraction::getAbsFrictionalStress() const
{
  pair<bool,double> res;

  // calculate distance
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  // check if there is contact
  if(dist<(eq_dist*eq_dist))
  {
    res.first=true;
    double Ac=eq_dist*eq_dist*0.7854; // contact area
    res.second=m_Ffric.norm()/Ac;
  }
  else
  {
    res.first=false;
  }

  return res;
}

/*!
  get max. frictional force, i.e. coeff. of friction * normal force
*/
pair<bool,double> CHertzMindlinViscoInteraction::getAbsMuFN() const
{
  pair<bool,double> res;

  // calculate distance
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  // check if there is contact
  if(dist<(eq_dist*eq_dist))
  { // contact -> calculate forces

    //--- viscoelastic force ---
    double R_ij=1.0/(1.0/m_p1->getRad()+1.0/m_p2->getRad());
    dist=sqrt(dist);
    double dn=eq_dist-dist;
    Vec3 dir=D.unit();
    double norm_force=4.0/3.0*m_E*sqrt(R_ij)*pow(dn,1.5);
    Vec3 force = dir*norm_force;

    double M_ij=1.0/(1.0/m_p1->getMass()+1.0/m_p2->getMass()); //equivalent mass
    Vec3 v_rel=m_p2->getVel()-m_p1->getVel(); //relative velocity
    Vec3 v_nor=(v_rel*dir)*dir; //normal relative velocity
    double lnCOR=log(m_COR);
    //see (Antypov,2011) eq.6:
    Vec3 normal_damping = 2.0*sqrt(5.0/3.0)*lnCOR/sqrt(pow(lnCOR,2)+pow(3.1415926,2))*sqrt(M_ij*m_E)*pow(R_ij,0.25)*pow(m_dn,0.25)*v_nor;
    force -= normal_damping;

    res.first=true;
    res.second=force.norm();
  }
  else
  {
    res.first=false;
  }

  return res;
}

/*!
  get max. frictional stress, i.e. coeff. of friction * normal stress
*/
pair<bool,double>
CHertzMindlinViscoInteraction::getMaxFricStress() const
{
  pair<bool,double> res;

  // calculate distance
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  // check if there is contact
  if(dist<(eq_dist*eq_dist))
  { // contact -> calculate forces
    //--- viscoelastic force ---
    double R_ij=1.0/(1.0/m_p1->getRad()+1.0/m_p2->getRad());
    dist=sqrt(dist);
    double dn=eq_dist-dist;
    Vec3 dir=D.unit();
    double norm_force=4.0/3.0*m_E*sqrt(R_ij)*pow(dn,1.5);
    Vec3 force = dir*norm_force;

    double M_ij=1.0/(1.0/m_p1->getMass()+1.0/m_p2->getMass()); //equivalent mass
    Vec3 v_rel=m_p2->getVel()-m_p1->getVel(); //relative velocity
    Vec3 v_nor=(v_rel*dir)*dir; //normal relative velocity
    double lnCOR=log(m_COR);
    //see (Antypov,2011) eq.6:
    Vec3 normal_damping = 2.0*sqrt(5.0/3.0)*lnCOR/sqrt(pow(lnCOR,2)+pow(3.1415926,2))*sqrt(M_ij*m_E)*pow(R_ij,0.25)*pow(m_dn,0.25)*v_nor;
    force -= normal_damping;

    res.first=true;
    double Ac=eq_dist*eq_dist*0.7854; // contact area
    res.second=force.norm()/Ac;
  }
  else
  {
    res.first=false;
  }

  return res;
}

/*!
	get current normal force
*/
pair<bool,double> CHertzMindlinViscoInteraction::getAbsFN() const
{
  return make_pair(m_is_touching,m_normal_force.norm());
}

/*!
  get current normal stress
*/
pair<bool,double>
CHertzMindlinViscoInteraction::getNormalStress() const
{
  pair<bool,double> res;

  if(m_is_touching){
    res.first=true;
    double eq_dist=m_p1->getRad()+m_p2->getRad();
    double Ac=eq_dist*eq_dist*0.7854; // contact area
    res.second=m_normal_force.norm()/Ac;
  } else {
    res.first=false;
  }
  return res;
}

/*!
  get "force deficit", i.e. the force needed to make the contact dynamic
*/

/*!
  get the slipping velocity, i.e. the absolute value of the
  tangential part of the relatve particle velocity
*/
pair<bool,double>
CHertzMindlinViscoInteraction::getSlipVelocity() const
{
  pair<bool,double> res;

  // calculate distance
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  if(D.norm()<=(m_p1->getRad()+m_p2->getRad()))
  { // if contact
    // normal
    Vec3 normal=D.unit();
    // relative velocity
    Vec3 v_rel=m_p2->getVel()-m_p1->getVel();
    // tangential component
    Vec3 v_tan=v_rel-(v_rel*normal)*normal;
    res=make_pair(true,v_tan.norm());
  } else {
    res.first=false;
  }
  return res;
}

/*!
  get the potential energy stored in the interaction
*/
double CHertzMindlinViscoInteraction::getPotentialEnergy() const
{
  const double e_pot_norm=0.5*m_normal_force*m_normal_force/m_E;
  return e_pot_norm;
}

/*!
  Get the static/dynamic status of the interaction. Returns 1 for a contact in
  dynamic friction, 0 for static or no contact
*/
double CHertzMindlinViscoInteraction::getSlipping() const
{
  const double res=m_is_slipping ? 1.0 : 0.0;
  return res;
}


/*!
  Get "sticking" contacts, i.e. return 1 if the contact is touching but not
  slipping, 0 otherwise
*/
double CHertzMindlinViscoInteraction::getSticking() const
{
  const double res=(m_is_touching && !m_is_slipping) ? 1.0 : 0.0;
  return res;
}

/*!
  return the amount of energy dissipated during the last time step
*/
double CHertzMindlinViscoInteraction::getDissipatedEnergy() const
{
  return m_E_diss;
}

/*!
  return the amount of normal viscous energy dissipated during the last time step
*/
double CHertzMindlinViscoInteraction::getNormalViscousEnergy() const
{
  return m_E_visc_normal;
}

/*!
  return the amount of shear viscous energy dissipated during the last time step
*/
double CHertzMindlinViscoInteraction::getShearViscousEnergy() const
{
  return m_E_visc_shear;
}

/*!
  get net force on particle1 imposed by this interaction.
  Returns Vec3::ZERO if particles are not in contact.
*/
Vec3 CHertzMindlinViscoInteraction::getForce() const
{
  const Vec3 f=m_is_touching ? m_Ffric-m_normal_force : Vec3(0.0,0.0,0.0);
  return f;
}

/*!
  If the particles are in contact, get normal force, if not in contact return
  (0,0,0)
*/
Vec3 CHertzMindlinViscoInteraction::getNormalForce() const
{
  const Vec3 f=m_is_touching ? m_normal_force : Vec3(0.0,0.0,0.0);
  return f;
}

/*!
  return 1 if particles are in contact, 0 otherwise
*/
double CHertzMindlinViscoInteraction::Count() const
{
  double res=m_is_touching ? 1.0 : 0.0;

  return res;
}

/*!
  Calculate effective coefficient of friction for this interaction for a given
  direction of the applied shear force. If the effective coefficient of
  friction is infinite, -1 is returned.

  \param dir the direction of the applied shear force
  \return the effective coefficient of friction if it is finite, -1 otherwise
    and -2 for no contact
*/
pair<bool,double> CHertzMindlinViscoInteraction::getMuEff(
  const Vec3& dir,
  const Vec3& norm
) const
{
  pair<bool,double> res;
  CParticle* p1;
  CParticle* p2;

  // sort particles, so that p1 is "above" the slip plane
  const Vec3 h=m_p1->getPos()-m_p2->getPos();
  if(h*norm>0.0)
  {
    p1=m_p1;
    p2=m_p2;
  }
  else
  {
    p1=m_p2;
    p2=m_p1;
  }
  // get contact normal
  Vec3 nc=p1->getPos()-p2->getPos();
  // get distance
  double dist=nc.norm();
  // check if contact
  if(dist<=(p1->getRad()+p2->getRad()))
  { // if contact
    // get direction of current movement
    Vec3 d=p1->getVel()-p2->getVel();
    // get tangential part
    d-=(d*nc.unit())*nc.unit();
    // calculate effective coefficient of friction
    double h1=(dir.unit()*d.unit())-m_mu*(dir.unit()*nc.unit());
    double h2=m_mu*(norm.unit()*nc.unit())+(norm.unit()*d.unit());
    if(h1>0)
    {
      res.first=true;
      res.second=h2/h1;
    }
    else
    {
      res.first=false;
    }
    cout << "positions : " << p1->getPos() << " , " << p2->getPos() << endl;
    cout << "velocities: " << p1->getVel() << " , " << p2->getVel() << endl;
    cout << "v_tan     : " << d << endl;
    cout << "h1,h2     : " << h1 << " , " << h2 << endl;
    cout << "mu_eff    : " << res.second << endl;
  }
  else
  {
    res.first=false;
  }
  return res;
}

/*!
  Get the particle member function which returns a scalar field of a given name.

  \param name the name of the field
*/
CHertzMindlinViscoInteraction::ScalarFieldFunction
CHertzMindlinViscoInteraction::getScalarFieldFunction(
  const string& name
)
{
  CHertzMindlinViscoInteraction::ScalarFieldFunction sf;

  if (name=="potential_energy")
  {
    sf=&CHertzMindlinViscoInteraction::getPotentialEnergy;
  }
  else if (name=="slipping")
  {
    sf=&CHertzMindlinViscoInteraction::getSlipping;
  }
  else if (name=="sticking")
  {
    sf=&CHertzMindlinViscoInteraction::getSticking;
  }
  else if (name=="count")
  {
    sf=&CHertzMindlinViscoInteraction::Count;
  }
  else if (name=="dissipated_energy")
  {
    sf=&CHertzMindlinViscoInteraction::getDissipatedEnergy;
  }
  else if (name=="viscous_energy_normal")
  {
    sf=&CHertzMindlinViscoInteraction::getNormalViscousEnergy;
  }
  else if (name=="viscous_energy_shear")
  {
    sf=&CHertzMindlinViscoInteraction::getShearViscousEnergy;
  }
  else
  {
    sf=NULL;
    cerr
      << "ERROR - invalid name for interaction scalar  access function"
      << endl;
  }

  return sf;
}

/*!
  Get the particle member function which returns a checked scalar field of a given name.

  \param name the name of the field
*/
CHertzMindlinViscoInteraction::CheckedScalarFieldFunction
CHertzMindlinViscoInteraction::getCheckedScalarFieldFunction(
  const string& name
)
{
  CHertzMindlinViscoInteraction::CheckedScalarFieldFunction sf;

  if (name=="mu_eff_xy")
  {
    sf=&CHertzMindlinViscoInteraction::getMuEffXY;
  }
  else if (name=="mu_eff_xz")
  {
    sf=&CHertzMindlinViscoInteraction::getMuEffXZ;
  }
  else if (name=="f_fric")
  {
    sf=&CHertzMindlinViscoInteraction::getAbsFrictionalForce;
  }
  else if (name=="fric_stress")
  {
    sf=&CHertzMindlinViscoInteraction::getAbsFrictionalStress;
  }
  else if (name=="f_normal")
  {
    sf=&CHertzMindlinViscoInteraction::getAbsFN;
  }
  else if (name=="normal_stress")
  {
    sf=&CHertzMindlinViscoInteraction::getNormalStress;
  }
  else if (name=="muF_n")
  {
    sf=&CHertzMindlinViscoInteraction::getAbsMuFN;
  }
  else if (name=="max_fric_stress")
  {
    sf=&CHertzMindlinViscoInteraction::getMaxFricStress;
  }
  else if (name=="v_slip")
  {
    sf=&CHertzMindlinViscoInteraction::getSlipVelocity;
  }
  else
  {
    sf=NULL;
    cerr
      << "ERROR - invalid name for interaction scalar  access function"
      << endl;
  }

  return sf;
}


/*!
  Get the particle member function which returns a vector field of a given name.

  \param name the name of the field
*/
CHertzMindlinViscoInteraction::VectorFieldFunction
CHertzMindlinViscoInteraction::getVectorFieldFunction(
  const string& name
)
{
  CHertzMindlinViscoInteraction::VectorFieldFunction vf;

  if (name=="force"){
    vf=&CHertzMindlinViscoInteraction::getForce;
  } else if (name=="normal_force") {
    vf = &CHertzMindlinViscoInteraction::getNormalForce;
  } else {
    vf=NULL;
    cerr
      << "ERROR - invalid name for interaction vector access function"
      << endl;
  }

  return vf;
}



/*!
  Pack a CHertzMindlinViscoInteraction into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface
::pack<CHertzMindlinViscoInteraction>(
  const CHertzMindlinViscoInteraction& I
)
{
  append(I.m_r0);
  append(I.m_E);
  append(I.m_nu);
  append(I.m_dn);
  append(I.m_mu);
  append(I.m_COR);
  append(I.m_dt);
  append(I.m_id[0]);
  append(I.m_id[1]);
  append(I.m_Ffric);
}


/*!
  Unpack a CHertzMindlinViscoInteraction from a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface
::unpack<CHertzMindlinViscoInteraction>(
  CHertzMindlinViscoInteraction& I
)
{
  I.m_r0=pop_double();
  I.m_E=pop_double();
  I.m_nu=pop_double();
  I.m_dn=pop_double();
  I.m_mu=pop_double();
  I.m_COR=pop_double();
  I.m_dt=pop_double();
  I.m_id.erase(I.m_id.begin(),I.m_id.end());
  I.m_id.push_back(pop_int());
  I.m_id.push_back(pop_int());
  I.m_Ffric=pop_vec3();
}


ostream& operator<<(
  ostream& ost,
  const CHertzMindlinViscoInteraction& BI
)
{
  ost << "[" << BI.m_p1->getID() << " - " << BI.m_p2->getID() << "]";
  return ost;
}
