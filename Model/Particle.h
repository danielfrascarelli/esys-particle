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

#ifndef __PARTICLE_H
#define __PARTICLE_H

// -- project includes --
#include "Foundation/Quaternion.h"
#include "Foundation/vec3.h"
#include "Foundation/Matrix3.h"
#include "Model/BasicParticle.h"
#include "Parallel/CheckPointable.h"

//--- STL includes ---
#include <map>
#include <vector>
#include <utility>
#include <string>
#include <iostream>

using std::map;
using std::vector;
using std::pair;
using std::string;

template <class T> class ParallelParticleArray;
class AMPISGBufferRoot;
class AMPIBuffer;

namespace esys
{
  namespace lsm
  {
    class SimpleParticleData;
  }
}

/*!
  \brief Class for a basic particle
*/
class CParticle : public CBasicParticle, public esys::lsm::CheckPointable
{
 public: // types
  class exchangeType
  {
  public:
    exchangeType()
      : m_pos(),
        m_initPos(),
        m_oldPos(),
      m_vel()
    {
      m_is_dynamic=true;
    }

      exchangeType(const Vec3 &pos, const Vec3 &initPos, const Vec3 &oldPos, const Vec3 &vel,bool is_dyn)
      : m_pos(pos),
        m_initPos(initPos),
        m_oldPos(oldPos),
        m_vel(vel)
    {
      m_is_dynamic=is_dyn;
    }
    
    Vec3 m_pos;
    Vec3 m_initPos;
    Vec3 m_oldPos;    
    Vec3 m_vel;
    bool m_is_dynamic;
  };
  
  typedef double (CParticle::* ScalarFieldFunction)() const; 
  typedef Vec3 (CParticle::* VectorFieldFunction)() const; 

 protected:
  //! stress tensor. \warning Warning: this is unscaled, i.e. without the 1/V term
  Matrix3 m_sigma; 
  Vec3 m_vel,m_force;
  Vec3 m_oldpos; //!< position at the time of last neighbor search
  Vec3 m_initpos; //!< position at time of construction
  Vec3 m_circular_shift; //!< shift vector if particle is circular image
  double m_mass,m_div_mass;
  double m_div_vol; //!< inverse of the particle volume, used in stress calculations

  bool flag;
  bool m_is_dynamic; 

  void setForce(const Vec3 &force) {m_force = force;}

 public:
  CParticle();
  CParticle(double,double,const Vec3&,const Vec3&,const Vec3&,int,bool);
  CParticle(double,double,const Vec3&,const Vec3&,const Vec3&,const Vec3&,const Vec3&,int,bool); // including oldpos
  CParticle(const esys::lsm::SimpleParticleData &particleData);
  virtual ~CParticle(){};


  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);
  
  inline const Vec3 &getInitPos() const {return m_initpos;}
  inline void setInitPos(const Vec3 &initPos) {m_initpos = initPos;}
  inline Vec3 getDisplacement() const {return (m_pos-m_oldpos);} ;
  inline Vec3 getTotalDisplacement() const {return (m_pos-m_initpos);} ;
  inline const Vec3 &getOldPos() const {return m_oldpos;};
  inline Vec3 getVel() const {return m_vel;};
  inline double getAbsVel() const {return m_vel.norm();};
  inline void setVel(const Vec3 &V){m_vel=V;};
  inline void setMass(double mass) {m_mass = mass; m_div_mass = 1.0/m_mass;}
  inline double getMass() const {return m_mass;};
  inline double getInvMass() const {return m_div_mass;};
  inline Vec3 getForce() const {return m_force;};
  virtual void setDensity(double); // needs to be virtual , different for rot. particle (mom. inert) 

  void resetDisplacement(){m_oldpos=m_pos;};
  double getIDField() const {return double(m_global_id);};
  double getTagField() const {return double(getTag());};
  void applyForce(const Vec3&,const Vec3&);
  virtual void integrate(double);
  virtual void integrateTherm(double){}
  virtual void zeroForce();
  virtual void zeroHeat() {}
  virtual void thermExpansion() {}
  inline void moveToRel(const Vec3 &v){m_pos=m_initpos+v;}; //!< move relative to initial position
  inline double getKineticEnergy() const {return 0.5*m_mass*m_vel*m_vel;};

  // switching on/off dynamic behaviour
  virtual void setNonDynamic() {m_is_dynamic=false;};
  virtual void setNonDynamicLinear() {m_is_dynamic=false;};
  virtual void setNonDynamicRot(){}; // do nothing
  virtual void resetRotation(){}; // do nothing
  inline void changeRadiusBy(double deltaR){m_rad += deltaR;}

  void setFlag(bool b=true){flag=b;};
  bool isFlagged() const {return flag;};
  void writeAsDXLine(ostream&,int slid=0);

  friend ostream& operator<<(ostream&, const CParticle&);
  void print(){cout << *this << endl << flush;};

  void rescale() {};
  exchangeType getExchangeValues();
  void setExchangeValues(const exchangeType&);

  // circular
  void setCircular(const Vec3&);

  // stress
  double sigma_xx_2D() const {return m_sigma(0,0)/(M_PI*m_rad*m_rad);};
  double sigma_xy_2D() const {return m_sigma(0,1)/(M_PI*m_rad*m_rad);};
  double sigma_yy_2D() const {return m_sigma(1,1)/(M_PI*m_rad*m_rad);};
  double sigma_d() const;
  double sigma_vonMises() const;
  Matrix3 sigma() const; 
	
  friend class TML_PackedMessageInterface;
  
  virtual void saveCheckPointData(std::ostream& oStream);
  virtual void saveSnapShotData(std::ostream& oStream);

  //virtual Quaternion getQuat(){return Quaternion(1.0,Vec3(0.0,0.0,0.0));};
  virtual void applyMoment(const Vec3&){};

  static void get_type() {cout <<" CParticle" ;};

  virtual void loadCheckPointData(std::istream &iStream);

  template <typename TmplVisitor>
  void visit(TmplVisitor &visitor)
  {
    visitor.visitParticle(*this);
  }

public:  
  // Ensure that particles only move in the x-y plane 2D computations
  inline static void setDo2dCalculations(bool do2dCalculations) {s_do2Calculations = do2dCalculations;}
  inline static bool getDo2dCalculations() {return s_do2Calculations;}

private:
  static bool s_do2Calculations;

 
};

/* CParticle extractCParticleFrom(AMPIBuffer*); */
/* CParticle extractCParticleFrom(AMPISGBufferRoot*,int); */

#endif //__PARTICLE_H

