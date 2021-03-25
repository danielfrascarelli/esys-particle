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

#include <math.h>

#include "Model/SpringDashpotFrictionInteraction.h"
#include "Foundation/console.h"
#include "tml/message/packed_message_interface.h"

CSpringDashpotFrictionIGP::CSpringDashpotFrictionIGP() : AIGParam(), youngModulus(0.0), mu(0.0), poissonRatio(0.0), cor(1.0), dt(0.0)
{
}

CSpringDashpotFrictionIGP::CSpringDashpotFrictionIGP(const std::string &name, double youngModulus, double poissionRatio, double cor, double fricCoef, double dT)
  : AIGParam(name),
    youngModulus(youngModulus),
    poissonRatio(poissonRatio),
    cor(cor),
    mu(fricCoef),
    dt(dT)
{
}

CSpringDashpotFrictionInteraction::CSpringDashpotFrictionInteraction()
{
  m_damping=0.0;
  m_k=0.0;
  m_mu=0.0;
  m_r0=0.0;
  m_ks=0.0;
  m_dt=0.0;
  m_is_slipping=false;
  m_is_touching=false;
  m_E_diss=0.0; 
}

void CSpringDashpotFrictionIGP::setTimeStepSize(double timeStepSize)
{
  this->dt = timeStepSize;
}

/*!
  constructor for CSpringDashpotFrictionInteraction without friction parameters, only calls
  the constructor of APairInteraction with the 2 particle pointers
*/
CSpringDashpotFrictionInteraction::CSpringDashpotFrictionInteraction(CParticle* p1,CParticle* p2):APairInteraction(p1,p2)
{
  m_is_slipping=false;
  m_is_touching=false;
  m_E_diss=0.0; 
}


CSpringDashpotFrictionInteraction::CSpringDashpotFrictionInteraction(CParticle* p1,CParticle* p2,const CSpringDashpotFrictionIGP& param):APairInteraction(p1,p2)
{
  // scale elastic param
  double f=1.0; 
  if(!CParticle::getDo2dCalculations()){
    f=0.25*M_PI*(p1->getRad()+p2->getRad());
  }

  m_k=f*param.youngModulus;
  m_mu=param.mu;
  
  m_damping=0.0;
  if ((param.cor < 1.0) && (param.cor >= 0.0)) {
    double A=log(param.cor);
    double M1 = p1->getMass();
    double M2 = p2->getMass();
    double Meff = (2.*M1*M2)/(M1+M2);
    m_damping = sqrt(m_k*Meff/(4./(A*A) - 1));
  }
  
  m_ks=f*param.youngModulus/2./(1.+param.poissonRatio);

  m_r0=p1->getRad()+p2->getRad();
  m_dt=param.dt;

  m_cpos=p1->getPos()+((p2->getPos()-p1->getPos())*p1->getRad()/m_r0);

  m_is_slipping=false;
  m_is_touching=false;
  m_E_diss=0.0; 
}

CSpringDashpotFrictionInteraction::~CSpringDashpotFrictionInteraction()
{
}

void CSpringDashpotFrictionInteraction::setTimeStepSize(double dt)
{
  m_dt = dt;
}

/*!
  Calculate elastic and frictional forces. 
*/
void CSpringDashpotFrictionInteraction::calcForces()
{
  // calculate distance
  const Vec3 D=m_p1->getPos()-m_p2->getPos();
  const double eq_dist=m_p1->getRad()+m_p2->getRad();
  double dist=D*D;
  // check if there is contact
  if(dist<(eq_dist*eq_dist))
  { // contact -> calculate forces
    //--- elastic force ---
    dist=sqrt(dist);
    double force_norm=m_k*(dist-eq_dist);

    Vec3 dV = m_p2->getVel() - m_p1->getVel();
    double relV = (dV*D)/dist;

    force_norm -= m_damping*relV;
    if (force_norm > 0.) {
       force_norm = 0.;
    }
    
    const Vec3 force=D*force_norm/dist; 
    
    m_normal_force=force;
    const Vec3 pos=m_p2->getPos()+(m_p2->getRad()/eq_dist)*D;
    // apply elastic force
    m_p2->applyForce(force,pos);
    m_p1->applyForce(-1.0*force,pos);
    //--- frictional force ---
    // particle movement since last timestep
    const Vec3 d1=m_p1->getVel()*m_dt;
    const Vec3 d2=m_p2->getVel()*m_dt;
    Vec3 ds = d2 - d1;
    // Compute tangential part by subtracting off normal component.
    ds -= ((ds*D)/(D.norm2()))*D;
    m_Ffric += m_ks * ds;
 
    const double FfricNorm = m_Ffric.norm();
    const double forceNorm = force.norm();
    // decide static/dynamic
    if (FfricNorm > forceNorm*m_mu)
    { // tangential force greater than static friction -> dynamic
      m_Ffric=m_Ffric*((m_mu*forceNorm)/FfricNorm);
      m_force_deficit=Vec3(0.0,0.0,0.0);
      m_is_slipping=true;
      m_E_diss = fabs(m_Ffric * ds); // energy dissipated
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

bool CSpringDashpotFrictionInteraction::isPersistent()
{
  const Vec3 D=m_p1->getPos()-m_p2->getPos();
  const double dist=D*D;
  const double eq_dist=m_p1->getRad()+m_p2->getRad();
  return (dist <= (eq_dist*eq_dist));
}

/*!
  get the force needed to overcome friction and make the interaction slip
*/
double CSpringDashpotFrictionInteraction::getAbsForceDeficit()const
{
  return m_force_deficit.norm();
}

/*!
  get current frictional/stopping force
*/
pair<bool,double> CSpringDashpotFrictionInteraction::getAbsFrictionalForce() const
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
pair<bool,double> CSpringDashpotFrictionInteraction::getAbsFrictionalStress() const
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
          double Ac=eq_dist*eq_dist*M_PI/4.0; // contact area
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
pair<bool,double> CSpringDashpotFrictionInteraction::getAbsMuFN() const
{
  pair<bool,double> res;

  // calculate distance
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  // check if there is contact
  if(dist<(eq_dist*eq_dist))
  {
    dist=sqrt(dist);
    Vec3 force=D*(m_k*(dist-eq_dist)/dist);
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
pair<bool,double> CSpringDashpotFrictionInteraction::getMaxFricStress() const
{
  pair<bool,double> res;

  // calculate distance
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  // check if there is contact
  if(dist<(eq_dist*eq_dist))
  {
    dist=sqrt(dist);
    Vec3 force=D*(m_k*(dist-eq_dist)/dist);
    res.first=true;
    double Ac=eq_dist*eq_dist*M_PI/4.0; // contact area
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
pair<bool,double> CSpringDashpotFrictionInteraction::getAbsFN() const
{
  return make_pair(m_is_touching,m_normal_force.norm());
}

/*!
	get current normal stress
*/
pair<bool,double> CSpringDashpotFrictionInteraction::getNormalStress() const
{
  pair<bool,double> res;

  if(m_is_touching){
    res.first=true;
    double eq_dist=m_p1->getRad()+m_p2->getRad();
    double Ac=eq_dist*eq_dist*M_PI/4.0; // contact area
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
pair<bool,double> CSpringDashpotFrictionInteraction::getSlipVelocity() const
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

  \warning For performance reasons the tangential part of the elastic energy is calculated directly from the current tangential force whereas it would be more accurate to calculate it incrementally during the force calculation. Be aware that therefore the potential energy is an approximation. Tests suggest an accuracy of a few percent.
*/
double CSpringDashpotFrictionInteraction::getPotentialEnergy() const
{
  const double e_pot_norm=0.5*m_normal_force*m_normal_force/m_k;
  double e_pot_tan=0.5*m_Ffric*m_Ffric/m_ks;
  
  return e_pot_norm+e_pot_tan;
}

/*!
  Get the static/dynamic status of the interaction. Returns 1 for a contact in dynamic
  friction, 0 for static or no contact
*/
double CSpringDashpotFrictionInteraction::getSlipping() const
{
  const double res=m_is_slipping ? 1.0 : 0.0;
  return res;
}


/*!
  Get "sticking" contacts, i.e. return 1 if the contact is touching but not 
  slipping, 0 otherwise
*/
double CSpringDashpotFrictionInteraction::getSticking() const
{
  const double res=(m_is_touching && !m_is_slipping) ? 1.0 : 0.0;
  return res;
}

/*!
  return the amount of energy dissipated during the last time step
*/
double CSpringDashpotFrictionInteraction::getDissipatedEnergy() const
{
  return m_E_diss;
}

/*!
  get net force on particle1 imposed by this interaction.
  Returns Vec3::ZERO if particles are not in contact.
*/
Vec3 CSpringDashpotFrictionInteraction::getForce() const
{
  const Vec3 f=m_is_touching ? m_Ffric-m_normal_force : Vec3(0.0,0.0,0.0);
  return f; 
}

/*!
  If the particles are in contact, get normal force, if not in contact return (0,0,0)
*/
Vec3 CSpringDashpotFrictionInteraction::getNormalForce() const
{
  const Vec3 f=m_is_touching ? m_normal_force : Vec3(0.0,0.0,0.0);
  return f; 
}

/*!
  return 1 if particles are in contact, 0 otherwise
*/
double CSpringDashpotFrictionInteraction::Count() const
{
  double res=m_is_touching ? 1.0 : 0.0;

  return res;
}

/*!
  Calculate effective coefficient of friction for this interaction for a given
  direction of the applied shear force. If the effective coefficient of friction
  is infinite, -1 is returned.

  \param dir the direction of the applied shear force
  \return the effective coefficient of friction if it is finite, -1 otherwise and -2 for no contact
*/
pair<bool,double> CSpringDashpotFrictionInteraction::getMuEff(const Vec3& dir,const Vec3& norm) const
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
CSpringDashpotFrictionInteraction::ScalarFieldFunction CSpringDashpotFrictionInteraction::getScalarFieldFunction(const string& name)
{
  CSpringDashpotFrictionInteraction::ScalarFieldFunction sf;

  if (name=="force_deficit")
  {
    sf=&CSpringDashpotFrictionInteraction::getAbsForceDeficit;
  }
  else if (name=="potential_energy")
  {
    sf=&CSpringDashpotFrictionInteraction::getPotentialEnergy;
  }
  else if (name=="slipping")
  {
    sf=&CSpringDashpotFrictionInteraction::getSlipping;
  }
  else if (name=="sticking")
  {
    sf=&CSpringDashpotFrictionInteraction::getSticking;
  }
  else if (name=="count")
  {
    sf=&CSpringDashpotFrictionInteraction::Count;
  }
  else if (name=="dissipated_energy")
  {
    sf=&CSpringDashpotFrictionInteraction::getDissipatedEnergy;
  }
  else
  {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction scalar  access function" << endl;
  }

  return sf;
}

/*!
  Get the particle member function which returns a checked scalar field of a given name.

  \param name the name of the field
*/
CSpringDashpotFrictionInteraction::CheckedScalarFieldFunction CSpringDashpotFrictionInteraction::getCheckedScalarFieldFunction(const string& name)
{
  CSpringDashpotFrictionInteraction::CheckedScalarFieldFunction sf;

  if (name=="mu_eff_xy")
  {
    sf=&CSpringDashpotFrictionInteraction::getMuEffXY;
  }
  else if (name=="mu_eff_xz")
  {
    sf=&CSpringDashpotFrictionInteraction::getMuEffXZ;
  }
  else if (name=="f_fric")
  {
    sf=&CSpringDashpotFrictionInteraction::getAbsFrictionalForce;
  }
  else if (name=="fric_stress")
  {
    sf=&CSpringDashpotFrictionInteraction::getAbsFrictionalStress;
  }
  else if (name=="f_normal")
  {
    sf=&CSpringDashpotFrictionInteraction::getAbsFN;
  }
  else if (name=="normal_stress")
  {
    sf=&CSpringDashpotFrictionInteraction::getNormalStress;
  }
  else if (name=="muF_n")
  {
    sf=&CSpringDashpotFrictionInteraction::getAbsMuFN;
  }
  else if (name=="max_fric_stress")
  {
    sf=&CSpringDashpotFrictionInteraction::getMaxFricStress;
  }
  else if (name=="v_slip")
  {
    sf=&CSpringDashpotFrictionInteraction::getSlipVelocity;
  }
  else
  {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction scalar  access function" << endl;
  }

  return sf;
}


/*!
  Get the particle member function which returns a vector field of a given name.
 
  \param name the name of the field
*/
CSpringDashpotFrictionInteraction::VectorFieldFunction CSpringDashpotFrictionInteraction::getVectorFieldFunction(const string& name)
{
  CSpringDashpotFrictionInteraction::VectorFieldFunction vf;

  if (name=="force"){
    vf=&CSpringDashpotFrictionInteraction::getForce;
  } else if (name=="normal_force") {
    vf = &CSpringDashpotFrictionInteraction::getNormalForce;    
  } else {
    vf=NULL;
    cerr << "ERROR - invalid name for interaction vector access function" << endl;
  }

  return vf;
}



/*!
  Pack a CSpringDashpotFrictionInteraction into a TML packed message
 
  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<CSpringDashpotFrictionInteraction>(const CSpringDashpotFrictionInteraction& I)
{
  append(I.m_k);
  append(I.m_r0);
  append(I.m_mu);
  append(I.m_ks);
  append(I.m_damping);
  append(I.m_dt);
  append(I.m_id[0]);
  append(I.m_id[1]);
  append(I.m_Ffric);
}

/*!
  Unpack a CSpringDashpotFrictionInteraction from a TML packed message
 
  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<CSpringDashpotFrictionInteraction>(CSpringDashpotFrictionInteraction& I)
{
  I.m_k=pop_double();
  I.m_r0=pop_double();
  I.m_mu=pop_double();
  I.m_ks=pop_double();
  I.m_damping=pop_double();
  I.m_dt=pop_double();
  I.m_id.erase(I.m_id.begin(),I.m_id.end());
  I.m_id.push_back(pop_int());
  I.m_id.push_back(pop_int());
  I.m_Ffric=pop_vec3();
}

/*!
  Save restart data to an open ostream 

  \param oStream the output stream
*/
void CSpringDashpotFrictionInteraction::saveRestartData(std::ostream &oStream)
{
  oStream << m_k << " ";
  oStream << m_r0 << " ";
  oStream << m_mu << " ";
  oStream << m_ks << " ";
  oStream << m_damping << " ";
  oStream << m_dt << " ";
  oStream << m_id[0] << " ";
  oStream << m_id[1] << " ";
  oStream << m_is_slipping << " ";
  oStream << m_is_touching << " ";
  oStream << m_Ffric.X() << " ";
  oStream << m_Ffric.Y() << " ";
  oStream << m_Ffric.Z() << std::endl;
}


/*!
  Load restart data from an open istream 

  \param iStream the input stream
*/
void CSpringDashpotFrictionInteraction::loadRestartData(std::istream &iStream)
{
  iStream >> m_k ;
  iStream >> m_r0 ;
  iStream >> m_mu ;
  iStream >> m_ks ;
  iStream >> m_damping ;
  iStream >> m_dt ;
  iStream >> m_id[0] ;
  iStream >> m_id[1] ;
  iStream >> m_is_slipping ;
  iStream >> m_is_touching ;
  iStream >> m_Ffric.X() ;
  iStream >> m_Ffric.Y() ;
  iStream >> m_Ffric.Z();
}

ostream& operator<<(ostream& ost,const CSpringDashpotFrictionInteraction& BI)
{
  ost << "[" << BI.m_p1->getID() << " - " << BI.m_p2->getID() << "]";
  return ost;
}
