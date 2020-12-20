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
#include "Foundation/console.h"
#include "Model/RotThermElasticInteraction.h"

CRotThermElasticIGP::CRotThermElasticIGP()
 : AIGParam(),
   m_kr(0.0),
   diffusivity(0.0)
{
}

CRotThermElasticIGP::CRotThermElasticIGP(
  const std::string &name,
  double normalK,
  double thermalDiffusivity
)
 : AIGParam(name),
   m_kr(normalK),
   diffusivity(thermalDiffusivity)
{
}

CRotThermElasticInteraction::CRotThermElasticInteraction() : ARotThermPairInteraction()
{
  m_kr          = 0.0;
  m_nForce      = 0.0;
  m_force       = Vec3(0.0,0.0,0.0);
  m_diffusivity = 0.0;
}

CRotThermElasticInteraction::CRotThermElasticInteraction(
  CRotThermParticle* p1,
  CRotThermParticle* p2,
  const CRotThermElasticIGP& param
)
 : ARotThermPairInteraction(p1,p2)
{
//  m_kr=param.m_kr;

// wyc added 22/02/2005
    double min_r ;
    double ran_ratio ;
    double ran_ratioH ;
   if (m_p1->getRad()<= m_p2->getRad() )    min_r = m_p1->getRad() ;
   else                                     min_r = m_p2->getRad() ;

//    double ran_ratio = 2.0*min_r/(m_p1->getRad()+m_p2->getRad());

   if(m_p1->getDo2dCalculations()) { // 2D
    ran_ratio  = 2.0*min_r/(m_p1->getRad()+m_p2->getRad());
    ran_ratioH = 2.0*min_r*(m_p1->getRad()+m_p2->getRad());
  }else{ // 3D
    ran_ratio  = 2.0*min_r*min_r/(m_p1->getRad()+m_p2->getRad());
    ran_ratioH = 2.0*min_r*min_r*(m_p1->getRad()+m_p2->getRad());
  }

  m_kr          = ran_ratio*param.m_kr;
  m_nForce      = 0.0;
  m_force       = Vec3(0.0,0.0,0.0);
  m_D           = p1->getPos()-p2->getPos();
  m_diffusivity = ran_ratioH*param.diffusivity;

}

Vec3 CRotThermElasticInteraction::getForce() const
{
/*
  Vec3 force = Vec3(0.0, 0.0, 0.0);
  Vec3 D = m_p1->getPos()-m_p2->getPos();
  double dist = D*D;
  double eq_dist = m_p1->getRad()+m_p2->getRad();
  if (dist < (eq_dist*eq_dist)) {
    dist = sqrt(dist);
    force = D*(m_kr*(dist - eq_dist)/dist);
  }
  return force;
*/
  return m_force;
}

/*!
  Calculate free elastic forces. 23 Flops if in contact, 10 Flops if not
 */
void CRotThermElasticInteraction::calcForces()
{ 
  //  console.XDebug() << "elastic interaction: [" << m_p1->getID() << " - " << m_p2->getID() << "]" << m_p1->getPos() << m_p2->getPos() << "\n"; 
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  //cout << "wyc in CRotEla ::calcF "  << m_p1->getID() << " - " << m_p2->getID()  <<endl;


  if(dist<(eq_dist*eq_dist)){ // contact -> calculate forces
    dist=sqrt(dist);
//    Vec3 force=D*m_kr*(dist-eq_dist)/dist;
    m_force=D.unit()*m_kr*(D.norm()-eq_dist);


//   cout << "dist= " << 2-D.norm() << " eq_dist= " <<2-eq_dist << endl;
//   cout << " D.unit=  " << D.unit() <<"   "<<D.norm()-eq_dist  << endl;  
//zjr 
//    cout << m_p1->getID() << " - " << m_p2->getID() << endl; 
//zjr    
//    cout << "wyc in elastic: " <<  force << endl ;
//zjr
//    cout << m_p1->getPos() << " - " << m_p2->getPos() << endl;


    Vec3 pos=m_p2->getPos()+(m_p2->getRad()/eq_dist)*D;
    m_p1->applyForce(-1.0*m_force,pos);
    m_p2->applyForce(m_force,pos);
    m_cpos=pos;
  } else { // no contact -> all forces are 0
    m_force = Vec3(0.0,0.0,0.0);
    m_nForce = 0.0;
  }
}

void CRotThermElasticInteraction::calcHeatTrans()
{

  double  eta = 1.5 ;

  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double Rij2 = D.norm2();
  double d_temp = m_p2->getTemperature() -  m_p1->getTemperature() ;
  double heatij = eta*m_diffusivity*d_temp/Rij2 ;

   m_p1->applyHeatTrans(heatij) ;
   m_p2->applyHeatTrans(-heatij) ;
}

/*!
  get the potential energy stored in the interaction
*/
double CRotThermElasticInteraction::getPotentialEnergy() const
{
//  const Vec3 normalForce = getForce();
//  return (0.5*normalForce*normalForce/m_kr);
  return (m_kr!=0.0) ? 0.5*(m_nForce*m_nForce)/m_kr : 0.0;
}

/*!
  Get the particle member function which returns a scalar field of a given name.
 
  \param name the name of the field
*/
CRotThermElasticInteraction::ScalarFieldFunction CRotThermElasticInteraction::getScalarFieldFunction(const string& name)
{
  CRotThermElasticInteraction::ScalarFieldFunction sf = NULL;
  
  if (name=="potential_energy")
  {
    sf=&CRotThermElasticInteraction::getPotentialEnergy;
  } else if (name=="count") {
    sf=&CRotThermElasticInteraction::Count;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction scalar access function" << endl;
  }
  
  return sf;
}

/*!
  Get the particle member function which returns a checked scalar field of a given name.

  \param name the name of the field
*/
CRotThermElasticInteraction::CheckedScalarFieldFunction CRotThermElasticInteraction::getCheckedScalarFieldFunction(const string& name)
{
  CRotThermElasticInteraction::CheckedScalarFieldFunction sf=NULL;

  cerr << "ERROR - invalid name for interaction vector access function" << endl;

  return sf;
}

/*!
  Get the particle member function which returns a vector field of a given name.
 
  \param name the name of the field
*/
CRotThermElasticInteraction::VectorFieldFunction CRotThermElasticInteraction::getVectorFieldFunction(const string& name)
{
  CRotThermElasticInteraction::VectorFieldFunction vf;
  
  if (name=="force"){
    vf=&CRotThermElasticInteraction::getForce;
  } else {
    vf=NULL;
    cerr << "ERROR - invalid name for interaction vector access function" << endl;
  }
  
  return vf;
}

/*!
  save restart data to ostream

  \param oStream the output stream
*/
void CRotThermElasticInteraction::saveRestartData(std::ostream &oStream)
{
  oStream << m_id[0] << " ";
  oStream << m_id[1] << " ";
  oStream << m_init  << " ";
  oStream << m_kr << " ";
  oStream << m_diffusivity << " ";
  oStream << m_D;
}

/*!
  load restart data from stream

  \param iStream the input stream
*/
void CRotThermElasticInteraction::loadRestartData(std::istream &iStream)
{
  iStream >> m_id[0];
  iStream >> m_id[1];
  iStream >> m_init ;
  iStream >> m_kr;
  iStream >> m_diffusivity;
  iStream >> m_D;
}

ostream& operator<<(ostream& ost,const CRotThermElasticInteraction& BI)
{
  ost << "[" << BI.m_p1->getID() << " - " << BI.m_p2->getID() << "]";
  return ost;
}
