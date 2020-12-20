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

#ifndef __THERMPARTICLE_H
#define __THERMPARTICLE_H

// -- project includes --
#include "Foundation/vec3.h"

// --- STL includes ---
#include <map>
#include <utility>


using std::map;
using std::pair;
using std::make_pair;

/*!
  Thermal Particle class.
*/
class CThermParticle
{
protected:

  double m_temperature ;
  double m_temperature_ini ;
  double m_Cp    ;
//  double m_density ;//   ???
  double m_heat_frict ;
  double m_heat_trans ;
  double m_therm_expansion0 ;
  double m_therm_expansion1 ;
  double m_therm_expansion2 ;
  double m_rad_ini ;    

public:
//  static const CBasicParticle INVALID;

  CThermParticle();
  CThermParticle(double rad_ini);
  CThermParticle(double temperature, 
                 double m_temperature_ini,
                 double Cp, 
                 double heat_frict, 
                 double heat_trans,
                 double therm_expansion0, 
                 double therm_expansion1, 
                 double therm_expansion2,
                 double rad_ini);
//  CThermParticle(const esys::lsm::SimpleParticleData &data);

  virtual ~CThermParticle(){};
/*
  inline Vec3 & getPPos() {return m_pos;}
  inline Vec3 getPos() const {return m_pos;}
  inline double getRad() const {return m_rad;}
  inline int getID() const {return m_global_id;}

  inline void moveBy(Vec3 v){m_pos+=v;} //!< move relative to current position
  inline void moveTo(Vec3 v){m_pos=v;} //!< move absolute
  inline void setRad(double r){m_rad=r;}

  //! particle tag handling
  inline void setTag(int t){m_tag=t;}
  inline int getTag() const {return m_tag;}
  inline bool isValid() const {return (getID() >= 0);}
};
*/

  inline void setTemperature(double t){m_temperature=t;}
  inline double getTemperature() const {return m_temperature;}

  inline void setEquilibTemperature(double t){m_temperature_ini=t;}
  inline double getEquilibTemperature() const {return m_temperature_ini;}

  inline void setEquilibRadius(double r){m_rad_ini=r;}
  inline double getEquilibRadius() const {return m_rad_ini;}

//  inline void setCp(double t){m_Cp = t;}
  inline double getCp() const {return m_Cp;}
  inline void setCp(double cp) {m_Cp = cp;}

  inline double getThermExpansion0() const {return m_therm_expansion0 ;}
  inline void setThermExpansion0(double te0) {m_therm_expansion0 = te0;}

  inline double getThermExpansion1() const {return m_therm_expansion1 ;}
  inline void setThermExpansion1(double te1) {m_therm_expansion1 = te1;}

  inline double getThermExpansion2() const {return m_therm_expansion2 ;}
  inline void setThermExpansion2(double te2) {m_therm_expansion2 = te2;}
//  void integrate_therm(double) ;
//  void zeroHeat() ;
//  void applyHeatTrans(const double) ;
  friend ostream& operator<<(ostream& ost,const CThermParticle& p);
};

ostream& operator<<(ostream&,const CThermParticle&);

#endif //__THERMPARTICLE_H
