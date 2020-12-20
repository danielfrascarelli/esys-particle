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

#include "ElasticInteraction.h"
#include "console.h"


CElasticIGP::CElasticIGP():AIGParam(), m_k(0.0)
{}

CElasticIGP::CElasticIGP(const std::string& name,double kn, bool scaling)
  : AIGParam(name), m_k(kn), m_scaling(scaling)
{}
  

CElasticInteraction::CElasticInteraction(
  CParticle* p1,
  CParticle* p2,
  const CElasticIGP& param
)
  : APairInteraction(p1,p2), m_is_touching(false)
{
  double effR=1.0;
  double effA=1.0; 
  double effL=1.0; // effective radius, cross section and length of the interaction for scaling 
  // equilibrium distance
  double r0=p1->getRad()+p2->getRad();
  
  m_scaling = param.m_scaling;
  // scale elastic param
  if (m_scaling) {
    if(!CParticle::getDo2dCalculations()){
      effR=0.5*r0;
    }
    effL=r0;
    effA = effR * effR;
  }
  m_k=param.m_k*effA/effL; 
  m_force=Vec3(0.0,0.0,0.0);
}

/*!
  Calculate free elastic forces. 23 Flops if in contact, 10 Flops if not
 */
void CElasticInteraction::calcForces()
{
  //  console.XDebug() << "elastic interaction: [" << m_p1->getID() << " - " << m_p2->getID() << "]" << m_p1->getPos() << m_p2->getPos() << "\n"; 
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  if(dist<(eq_dist*eq_dist)){
    dist=sqrt(dist);
    m_force=D*(m_k*(dist-eq_dist)/dist);
    Vec3 pos=m_p2->getPos()+(m_p2->getRad()/eq_dist)*D;
    m_p2->applyForce(m_force,pos);
    m_p1->applyForce(-1.0*m_force,pos); 
    m_cpos=pos;
    m_is_touching=true;
  } else {
    m_is_touching=false;
  }
}

/*!
  get the potential energy stored in the interaction
*/
double CElasticInteraction::getPotentialEnergy() const
{
  double e_pot_norm=0.5*m_force*m_force/m_k;

  return e_pot_norm;
}

Vec3 CElasticInteraction::getForce() const
{
  return m_force;
}

/*!
  return 1 if particles are in contact, 0 otherwise
*/
double CElasticInteraction::Count() const
{
  double res=m_is_touching ? 1.0 : 0.0;

  return res;
}

void CElasticInteraction::setStiffness(double k)
{
	m_k=k;
	console.XDebug() << "CElasticInteraction stiffness set to " << k << "\n";
}

/*!
  Get the particle member function which returns a scalar field of a given name.
 
  \param name the name of the field
*/
CElasticInteraction::ScalarFieldFunction CElasticInteraction::getScalarFieldFunction(const string& name)
{
  CElasticInteraction::ScalarFieldFunction sf;

  if (name=="potential_energy")
    {
      sf=&CElasticInteraction::getPotentialEnergy;
    }
  else if (name=="count")
    {
      sf=&CElasticInteraction::Count;
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
CElasticInteraction::CheckedScalarFieldFunction CElasticInteraction::getCheckedScalarFieldFunction(const string& name)
{
  CElasticInteraction::CheckedScalarFieldFunction sf;

  sf=NULL;
  cerr << "ERROR - invalid name for interaction scalar  access function" << endl;
  
  return sf;
}

/*!
  Get the particle member function which returns a vector field of a given name.
 
  \param name the name of the field
*/
CElasticInteraction::VectorFieldFunction CElasticInteraction::getVectorFieldFunction(const string& name)
{
  CElasticInteraction::VectorFieldFunction vf;

  if (name=="force"){
    vf=&CElasticInteraction::getForce;
  } else {
    vf=NULL;
    cerr << "ERROR - invalid name for interaction vector access function" << endl;
  }
  
  return vf;
}

/*!
	Get the interaction member function which sets a parameter of a given name
	
	\param pname the parameter name
*/
CElasticInteraction::ScalarSetFunction CElasticInteraction::getScalarSetFunction(const string& pname)
{
	CElasticInteraction::ScalarSetFunction sf=NULL;
	
	if (pname=="k") {
		sf=&CElasticInteraction::setStiffness;
	} else {
		console.Error() << "invalid name for setter function\n";
	}
	
	return sf; 
}

/*!
  save restart data to ostream

  \param oStream the output stream
*/
void CElasticInteraction::saveRestartData(std::ostream &oStream)
{
  oStream << m_id[0] << " ";
  oStream << m_id[1] << " ";
  oStream << m_init  << " ";
  oStream << m_k << " ";
  oStream << m_scaling;
}

/*!
  load restart data from stream

  \param iStream the input stream
*/
void CElasticInteraction::loadRestartData(std::istream &iStream)
{
  iStream >> m_id[0];
  iStream >> m_id[1];
  iStream >> m_init ;
  iStream >> m_k;
  iStream >> m_scaling;
}

ostream& operator<<(ostream& ost,const CElasticInteraction& BI)
{
  ost << "[" << BI.m_p1->getID() << " - " << BI.m_p2->getID() << "]";
  return ost;
}

