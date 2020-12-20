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

#include "VWFrictionInteraction.h"
#include "tml/message/packed_message_interface.h"

VWFrictionIGP::VWFrictionIGP() : CFrictionIGP(), m_alpha(0.0)
{}

VWFrictionIGP::VWFrictionIGP(const std::string &name, double kn, double mu, double ks, double dt, double a) 
  : CFrictionIGP(name,kn,mu,ks,dt)
{
  m_alpha=a;
}

/*!
	Default constructor
*/
CVWFriction::CVWFriction():CFrictionInteraction()
{
  m_alpha=0.0;
}

/*!
*/
CVWFriction::CVWFriction(CParticle* p1,CParticle* p2,const VWFrictionIGP& Param):CFrictionInteraction(p1,p2,Param)
{
	m_alpha=Param.m_alpha;
}

CVWFriction::~CVWFriction()
{}

/*!
	Calculate elastic & frictional forces.
	The current coefficient of friction is calculated by a velocity weakening
	friction law \f$ \mu=\frac{\mu_0}{1+2\alpha |v|} \f$
*/
void CVWFriction::calcForces()
{
	Vec3 pos;
	Vec3 force;
	double mu; // current coefficent of friction

	// calculate distance
	Vec3 D=m_p1->getPos()-m_p2->getPos();
	double dist=D*D;
	double eq_dist=m_p1->getRad()+m_p2->getRad();
	// check if there is contact
	if(dist<(eq_dist*eq_dist))
	{ // contact -> calculate forces
		//--- elastic force ---
		dist=sqrt(dist);
		force=D*(m_k*(dist-eq_dist)/dist);
		m_normal_force=force;
		pos=m_p2->getPos()+(m_p2->getRad()/eq_dist)*D;
		// apply elastic force
		m_p2->applyForce(force,pos);
		m_p1->applyForce(-1.0*force,pos);
		//--- frictional force ---
		// particle movement since last timestep
		Vec3 d1=m_p1->getVel()*m_dt;
		Vec3 d2=m_p2->getVel()*m_dt;
		Vec3 dFF=m_ks*(d2-d1);
		// tangential part
		Vec3 normal=D.unit();
		dFF-=(dFF*normal)*normal;
		m_Ffric+=dFF;
		// --- calc current mu ---
		// relative velocity
		Vec3 v_rel=m_p2->getVel()-m_p1->getVel();
		// tangential component
		Vec3 v_tan=v_rel-(v_rel*normal)*normal;
		// mu
		mu=m_mu/(1.0+2.0*m_alpha*v_tan.norm());
		// decide static/dynamic
		if(m_Ffric.norm()>force.norm()*mu)
		{ // tangential force greater than static friction -> dynamic
			m_Ffric=m_Ffric*((mu*force.norm())/m_Ffric.norm());
			m_force_deficit=Vec3(0.0,0.0,0.0);
			m_is_slipping=true;
		}
		else if(m_Ffric.norm()>0.0)
		{ // static friction
			m_is_slipping=false;
		}
		else
		{ // no frictional force -> force deficit=mu*F_n
			m_is_slipping=false;
		}
		m_p1->applyForce(m_Ffric,pos);
		m_p2->applyForce(-1.0*m_Ffric,pos);
		m_cpos=pos;
		m_is_touching=true;
	}
	else
	{ // no contact -> all forces are 0
		m_Ffric=Vec3(0.0,0.0,0.0);
		m_normal_force=Vec3(0.0,0.0,0.0);
		m_is_slipping=false;
		m_is_touching=false;
	}
}

/*!
	get current coefficient of friction
	*/
pair<bool,double> CVWFriction::getCurrentMu() const
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
		// mu
		double mu=m_mu/(1.0+2.0*m_alpha*v_tan.norm());
		res=make_pair(true,mu);
	} else {
		res.first=false;
	}
	return res;
}


CVWFriction::ScalarFieldFunction CVWFriction::getScalarFieldFunction(const string& name)
{
  CVWFriction::ScalarFieldFunction sf;

  if (name=="potential_energy"){
    sf=&CVWFriction::getPotentialEnergy;
  } else if (name=="slipping"){
    sf=&CVWFriction::getSlipping;
  } else if (name=="sticking"){
    sf=&CFrictionInteraction::getSticking;
  } else if (name=="count"){
    sf=&CVWFriction::Count;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction scalar access function" << endl;
  }

  return sf;
}

CVWFriction::CheckedScalarFieldFunction CVWFriction::getCheckedScalarFieldFunction(const string& name)
{
	CVWFriction::CheckedScalarFieldFunction sf;

	if (name=="mu_eff_xy"){
		sf=&CVWFriction::getMuEffXY;
	} else if (name=="mu_eff_xz") {
		sf=&CVWFriction::getMuEffXZ;
	} else if (name=="f_fric") {
		sf=&CVWFriction::getAbsFrictionalForce;
	} else if (name=="muF_n") {
		sf=&CVWFriction::getAbsMuFN;
	} if (name=="mu_current") {
		sf=&CVWFriction::getCurrentMu;
	} else if (name=="v_slip") {
		sf=&CVWFriction::getSlipVelocity;
	} else {
		sf=NULL;
		cerr << "ERROR - invalid name for checked interaction scalar access function" << endl;
	}

	return sf;
	}

CVWFriction::VectorFieldFunction CVWFriction::getVectorFieldFunction(const string&)
{
  CVWFriction::VectorFieldFunction vf;

	vf=NULL;

	return vf;
}


/*!
  Pack a CFrictionInteraction into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<CVWFriction>(const CVWFriction& I)
{
	append(I.m_k);
	append(I.m_r0);
	append(I.m_mu);
	append(I.m_ks);
	append(I.m_dt);
	append(I.m_alpha);
	append(I.m_id[0]);
	append(I.m_id[1]);
}

/*!
  Unpack a CFrictionInteraction from a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<CVWFriction>(CVWFriction& I)
{
	I.m_k=pop_double();
	I.m_r0=pop_double();
	I.m_mu=pop_double();
	I.m_ks=pop_double();
	I.m_dt=pop_double();
	I.m_alpha=pop_double();
	I.m_id.erase(I.m_id.begin(),I.m_id.end());
	I.m_id.push_back(pop_int());
	I.m_id.push_back(pop_int());
}

