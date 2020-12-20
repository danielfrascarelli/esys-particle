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

#ifndef __ROTPARTICLEVI_H
#define __ROTPARTICLEVI_H

// -- project includes --
#include "Foundation/vec3.h"
#include "Foundation/Matrix3.h"
#include "Model/Particle.h"
#include "Foundation/Quaternion.h"

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
  Class for a rotational particle, Verlet integration
*/
class CRotParticleVi : public CParticle
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
        m_quat()
    {
    }

    exchangeType(
      const Vec3 &pos,
      const Vec3 &initPos,
      const Vec3 &vel,
      const Vec3 &AngVel,
      const Vec3 &currAngVel,
      const Quaternion &quat
    )
      : m_pos(pos),
        m_initPos(initPos),
        m_vel(vel),
        m_angVel(AngVel),
        m_angVel_t(currAngVel),
        m_quat(quat)
    {
    }
  public:
    Vec3       m_pos;
    Vec3       m_initPos;
    Vec3       m_vel;
    Vec3       m_angVel;
    Vec3       m_angVel_t ;
    Quaternion m_quat;

    friend class TML_PackedMessageInterface;    
  };
  typedef double (CRotParticleVi::* ScalarFieldFunction)() const;
  typedef Vec3 (CRotParticleVi::* VectorFieldFunction)() const;

 protected:
 
  Quaternion  m_quat;
  Quaternion  m_initquat;
  Vec3        m_angVel;    // !  Angular velocity at time t -0.5*dt
  Vec3        m_angVel_t ; //! Angular velocity at time t
  Vec3        m_moment;
  double      m_inertRot;
  double      m_div_inertRot;
  bool        m_is_dynamic; 

 public:
  CRotParticleVi();
  CRotParticleVi(const esys::lsm::SimpleParticleData &particleData);
  CRotParticleVi(
    double      rad,
    double      mass,
    const Vec3& pos,
    const Vec3& vel,
    const Vec3& force,
    int         id,
    bool        is_dyn
  );
  CRotParticleVi(
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
    const Vec3& angvel_t  
  );
  CRotParticleVi(
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
    const Vec3&       angvel_t  
  );
  
  CRotParticleVi(const CParticle &p);

  virtual ~CRotParticleVi(){};

  static int getPackSize();
  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  // Need to define this for template class using forAllParticles call in Parallel/SubLattice.hpp.
  Vec3 getDisplacement() const {return CParticle::getDisplacement();};
  void resetDisplacement() {CParticle::resetDisplacement();}

  inline const Vec3 &getAngVel () const { return m_angVel;}
  inline const Vec3 getAngVel_t () const { return m_angVel_t ;} ;
  inline void setAngVel_t (const Vec3 &v) { m_angVel_t = v;} ;
  inline Vec3 getAngVelNR () const { return m_angVel_t;} // for field functions 
  inline void setAngVel(const Vec3 &V) { m_angVel = V;}  
  inline Quaternion getInitQuat() const { return m_initquat;}
  inline Quaternion getQuat() const { return m_quat;}
  inline void setQuat(const Quaternion &q) {m_quat = q;}
  inline double getInertRot () const { return m_inertRot; }
  inline void setInertRot (double inertRot)
  {
    m_inertRot = inertRot;
    m_div_inertRot = 1.0/m_inertRot;
  }
  inline double getInvInertRot () const { return m_div_inertRot; }
  inline Vec3 getMoment()  const {return m_moment;}
  inline void setMoment(const Vec3 &moment)  {m_moment = moment;}
  Vec3 getAngVector()const ;
  void applyMoment( const Vec3&);
  void integrate(double);
  void zeroForce();

//wycnewadded
  virtual void zeroHeat(){} ;
  virtual void integrateTherm(double){} ;
  virtual void setTemperature(double){} ;
  virtual void setCp(double){} ;
  virtual void setThermExpansion0(double){} ;
  virtual void setThermExpansion1(double){} ;
  virtual void setThermExpansion2(double){} ;
  virtual void thermExpansion(){} ;
  virtual double get_y() {return m_pos.Y(); } ;


  void rescale();
  void setCircular(const Vec3& cv);
  double getAngularKineticEnergy() const {return (0.2*m_mass*m_rad*m_rad)*(m_angVel_t*m_angVel_t);}
  double getLinearKineticEnergy() const {return (0.5*m_mass)*(m_vel*m_vel);}
  double getKineticEnergy() const {return getLinearKineticEnergy() + getAngularKineticEnergy();}
  void writeAsDXLine(ostream&,int slid=0);
  virtual void setNonRot() {m_inertRot=0.0;m_div_inertRot=0.0;};

  inline Quaternion getQuatFromRotVec(const Vec3 &vec) const
  {
    const double angle     = vec.norm();
    const double halfAngle = angle/2.0;
    if (halfAngle > 0.0) {
      return Quaternion(cos(halfAngle), (vec)*(sin(halfAngle)/angle));
    }
    return Quaternion();
  }
  void rotateBy(const Vec3 &vec) {m_quat = getQuatFromRotVec(vec)*m_quat;}
  void rotateTo(const Vec3 &vec) {m_quat = getQuatFromRotVec(vec);}

  static map<string,AField*> generateFields(ParallelParticleArray<CRotParticleVi>*);

  friend ostream& operator<<(ostream&, const CRotParticleVi&);
  void print(){cout << *this << endl << flush;};

  // -- checkpointing --
  virtual void saveSnapShotData(std::ostream& oStream);
  virtual void saveCheckPointData(std::ostream& oStream);
  virtual void loadCheckPointData(std::istream& iStream);

  CRotParticleVi::exchangeType getExchangeValues();
  void setExchangeValues(const CRotParticleVi::exchangeType &e);

  // stress 
  inline double sigma_xx_2D() const {return m_sigma(0,0)/(M_PI*m_rad*m_rad);};
  inline double sigma_xy_2D() const {return m_sigma(0,1)/(M_PI*m_rad*m_rad);};
  inline double sigma_yy_2D() const {return m_sigma(1,1)/(M_PI*m_rad*m_rad);};
//  inline double sigma_d() const;
  static void get_type() {cout <<" CRotParticleVi" ;};
  friend class TML_PackedMessageInterface;

  template <typename TmplVisitor>
  void visit(TmplVisitor &visitor)
  {
    visitor.visitRotParticleVi(*this);
  }

};

#endif //__ROTPARTICLEVI_H
