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

#include "Model/ThermParticle.h"

CThermParticle::CThermParticle()
  : m_temperature(0.0),
     m_temperature_ini(0.0),
    m_Cp(0.0),
    m_heat_frict(0.0),
    m_heat_trans(0.0),
    m_therm_expansion0(0.0),
    m_therm_expansion1(0.0),
    m_therm_expansion2(0.0),
    m_rad_ini(0.0)
{
}

CThermParticle::CThermParticle(double rad_ini)
  : m_temperature(0.0),
    m_temperature_ini(0.0),
    m_Cp(0.0),
    m_heat_frict(0.0),
    m_heat_trans(0.0),
    m_therm_expansion0(0.0),
    m_therm_expansion1(0.0),
    m_therm_expansion2(0.0),
    m_rad_ini(rad_ini)
{
}

CThermParticle::CThermParticle(double temperature, 
                               double temperature_ini,
                               double Cp, 
                               double heat_frict, 
                               double heat_trans,
                               double therm_expansion0,
                               double therm_expansion1,
                               double therm_expansion2,
                               double rad_ini)
  : m_temperature(temperature),
    m_temperature_ini(temperature_ini),
    m_Cp(Cp),
    m_heat_frict(heat_frict),
    m_heat_trans(heat_trans),
    m_therm_expansion0(therm_expansion0),
    m_therm_expansion1(therm_expansion1),
    m_therm_expansion2(therm_expansion2),
    m_rad_ini(rad_ini)
{
}

ostream& operator<<(ostream& ost,const CThermParticle& p)
{
  ost
    << "m_temperature=" << p.m_temperature
    << ", m_temperature_ini=" << p.m_temperature_ini
    << ", m_Cp" << p.m_Cp
    << ", m_heat_frict=" << p.m_heat_frict
    << ", m_heat_trans=" << p.m_heat_trans
    << ", m_therm_expansion0=" << p.m_therm_expansion0
    << ", m_therm_expansion1=" << p.m_therm_expansion1
    << ", m_therm_expansion2=" << p.m_therm_expansion2
    << ", m_rad_ini=" << p.m_rad_ini;

  return ost;
}
