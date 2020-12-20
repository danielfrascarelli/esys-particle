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

#ifndef __ROTTHERMPARTICLE_H
#define __ROTTHERMPARTICLE_H

// -- project includes --
#include "Foundation/vec3.h"
#include "Foundation/Matrix3.h"
#include "Model/ThermParticle.h"
#include "Foundation/Quaternion.h"
#include "Model/RotParticleVi.h"

template <class T> class ParallelParticleArray;
class AMPISGBufferRoot;
class AMPIBuffer;
class AField;

//--- MPIincludes ---
#include <mpi.h>

//--- STL includes ---
#include <map>
#include <vector>
#include <utility>
#include <string>

using std::map;
using std::vector;
using std::pair;
using std::string;

namespace esys
{
  namespace lsm
  {
    class SimpleParticleData;
  }
}


/*!
  Thermal Particle class.
*/
class CRotThermParticle : public CRotParticleVi, public CThermParticle
{

 public: // types

  class exchangeType
  {
  public:
    exchangeType()
      : m_pos(),
        m_initPos(),
        m_vel(),
        m_angVel(),
        m_angVel_t(),
        m_quat(),
        m_temperature(),
        m_temperature_ini() 
    {
    }

    exchangeType(
      const Vec3 &pos,
      const Vec3 &initPos,
      const Vec3 &vel,
      const Vec3 &AngVel,
      const Vec3 &currAngVel,
      const Quaternion &quat,
      const double temperature,
      const double temperature_ini
    )
      : m_pos(pos),
        m_initPos(initPos),
        m_vel(vel),
        m_angVel(AngVel),
        m_angVel_t(currAngVel),
        m_quat(quat),
        m_temperature(temperature),
        m_temperature_ini(temperature_ini)
    {
    }
  public:
    Vec3       m_pos;
    Vec3       m_initPos;
    Vec3       m_vel;
    Vec3       m_angVel;
    Vec3       m_angVel_t ;
    Quaternion m_quat;
    double m_temperature ;
    double m_temperature_ini ;

    friend class TML_PackedMessageInterface;
  };
  typedef double (CRotThermParticle::* ScalarFieldFunction)() const;
  typedef Vec3 (CRotThermParticle::* VectorFieldFunction)() const;





protected:

//  double m_tempa ;
//  double m_Cp    ;
//  double m_density ;
//  double m_heat_frict ;
//  double m_heat_trans ;
  

public:
//  static const CBasicParticle INVALID;

  CRotThermParticle();
  CRotThermParticle(const esys::lsm::SimpleParticleData &data);

  CRotThermParticle(const CRotParticleVi &p);

  CRotThermParticle(const CParticle &p);

  CRotThermParticle(
    double      rad,
    double      mass,
    const Vec3& pos,
    const Vec3& vel,
    const Vec3& force,
    int         id,
    bool        is_dyn
  );

  CRotThermParticle(
    double      rad,
    double      mass,
    const Vec3& pos,
    const Vec3& vel,
    const Vec3& force,
    int         id,
    Quaternion& quat,
    double      inertRot,
    const Vec3& moment,
    const Vec3& angvel,
    const Vec3& angvel_t,
    double temperature,
    double temperature_ini,
    double Cp,
    double heat_frict,
    double heat_trans,
    double therm_expansion0,
    double therm_expansion1,
    double therm_expansion2 
  );
  CRotThermParticle(
    double            rad,
    double            mass,
    const Vec3&       pos,
    const Vec3&       oldpos,
    const Vec3&       initpos,
    const Vec3&       vel,
    const Vec3&       force,
    int               id,
    const Quaternion& quat,
    const Quaternion& initquat,
    double            inertRot,
    const Vec3&       moment,
    const Vec3&       angvel,
    const Vec3&       angvel_t,
    double temperature,
     double temperature_ini,
    double Cp,
    double heat_frict,
    double heat_trans,
    double therm_expansion0,
    double therm_expansion1,
    double therm_expansion2

  );


   ~CRotThermParticle(){};

  void applyHeatTrans(const double); 
  void applyHeatFrict(const double);
  void integrateTherm(double);
  void zeroHeat();
  void thermExpansion() ;

  void integrate(double);
  inline void setTemperature(double t){m_temperature=t; m_temperature_ini=t;} ;
//  inline void setTemperatureIni(double t){m_temperature_ini=t;} ;
  inline double get_y() {return m_pos.Y(); } ;
  inline void setCp(double t) {m_Cp=t; } ; 
  inline void setThermExpansion0(double t) { m_therm_expansion0=t;} ;
  inline void setThermExpansion1(double t) { m_therm_expansion1=t;} ;
  inline void setThermExpansion2(double t) { m_therm_expansion2=t;} ;  

  void setCircular(const Vec3& cv);

  Vec3 getDisplacement() const {return CParticle::getDisplacement();};
  void resetDisplacement() {CParticle::resetDisplacement();}



  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  static map<string,AField*> generateFields(ParallelParticleArray<CRotThermParticle>*);

  friend ostream& operator<<(ostream&, const CRotThermParticle&);
  void print(){cout << *this << endl << flush;};

//  virtual void saveCheckPointData(std::ostream& oStream);
//  virtual void loadCheckPointData(std::istream& iStream);

  CRotThermParticle::exchangeType getExchangeValues();
  void setExchangeValues(const CRotThermParticle::exchangeType &e);


  static void get_type() {cout <<" CRotThermParticle" ;};
  friend class TML_PackedMessageInterface;

  template <typename TmplVisitor>
  void visit(TmplVisitor &visitor)
  {
    visitor.visitRotThermParticle(*this);
  }

};

// ostream& operator<<(ostream&,const CRotThermParticle&);

#endif //__ROTTHERMPARTICLE_H
