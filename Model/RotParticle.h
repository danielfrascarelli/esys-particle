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

#ifndef __ROTPARTICLE_H
#define __ROTPARTICLE_H

//--- MPIincludes ---
#include <mpi.h>

// -- project includes --
#include "Foundation/vec3.h"
#include "Foundation/Matrix3.h"
#include "Model/Particle.h"
#include "Foundation/Quaternion.h"

template <class T> class ParallelParticleArray;
class AMPISGBufferRoot;
class AMPIBuffer;


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
  \brief Class for a rotational particle

*/
class CRotParticle : public CParticle
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
        m_quat(),
        m_stress()
    {
      m_is_dynamic=true;
      m_is_rot=true;
    }

    exchangeType(
      const Vec3 &pos,
      const Vec3 &initPos,
      const Vec3 &vel,
      const Vec3 &currAngVel,
      const Quaternion &quat,
      const bool is_dyn,
      const bool is_rot,
      const Matrix3 stress
    )
      : m_pos(pos),
        m_initPos(initPos),
        m_vel(vel),
        m_angVel(currAngVel),
        m_quat(quat),
        m_stress(stress)
    {
      m_is_dynamic=is_dyn;
      m_is_rot=is_rot;
    }
  public:
    Vec3       m_pos;
    Vec3       m_initPos;
    Vec3       m_vel;
    Vec3       m_angVel;
    Quaternion m_quat;
    bool m_is_dynamic;
    bool m_is_rot;
    Matrix3 m_stress;
  
    friend class TML_PackedMessageInterface;    
  };
  typedef double (CRotParticle::* ScalarFieldFunction)() const;
  typedef Vec3 (CRotParticle::* VectorFieldFunction)() const;

 protected:
 
  Quaternion  m_quat;
  Quaternion  m_initquat;
  Vec3        m_angVel; //! Angular velocity at time t
  Vec3        m_moment;
  double      m_inertRot;
  double      m_div_inertRot;
  bool m_is_rot; //! false if rotational dynamics are switched off 

  void setMoment(const Vec3 &moment) {m_moment = moment;}

 public:
  CRotParticle();
  CRotParticle(const esys::lsm::SimpleParticleData &particleData);
  CRotParticle(const CParticle &particle);
  CRotParticle(
    double      rad,
    double      mass,
    const Vec3& pos,
    const Vec3& vel,
    const Vec3& force,
    int         id,
    bool is_dyn,
    bool is_rot
  );
  CRotParticle(
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
    bool is_dyn, 
    bool is_rot
  );
  CRotParticle(
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
    bool is_dyn,
    bool is_rot
  );
  
  virtual ~CRotParticle(){};

  static int getPackSize();
  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  // Need to define this for template class using forAllParticles call in Parallel/SubLattice.hpp.
  Vec3 getDisplacement() const {return CParticle::getDisplacement();};
  void resetDisplacement() {CParticle::resetDisplacement();}
  virtual void setDensity(double);
  
  inline const Vec3 &getAngVel () const { return m_angVel;}
  inline Vec3 getAngVelNR () const { return m_angVel;} // for field functions 
  inline void setAngVel(const Vec3 &V) { m_angVel = V;}  
  inline Quaternion getInitQuat() const { return m_initquat;}
  inline Quaternion getQuat() const { return m_quat;}
  inline void setQuat(const Quaternion &quat) { m_quat = quat;}
  inline double getInertRot () const { return m_inertRot; }
  inline void setInertRot (double inertRot)
  {
    m_inertRot = inertRot;
    m_div_inertRot = 1.0/m_inertRot;
  }
  inline double getInvInertRot () const { return m_div_inertRot; }
  inline Vec3 getMoment()  const {return m_moment;}
  void applyMoment( const Vec3&);
  void integrate(double);
  void integrateTherm(double dt){}
  virtual void thermExpansion() {}
  void zeroForce();
  virtual void zeroHeat(){}
  void rescale();
  void setCircular(const Vec3& cv);
  double getAngularKineticEnergy() const {return 0.5*m_inertRot*m_angVel*m_angVel;} // 1/2 I*omega^2
  double getLinearKineticEnergy() const {return (0.5*m_mass)*(m_vel*m_vel);}
  double getKineticEnergy() const {return getLinearKineticEnergy() + getAngularKineticEnergy();}
  void writeAsDXLine(ostream&,int slid=0);
  
  // switching on/off dynamic behaviour
  virtual void setNonDynamic() {m_is_dynamic=false;m_is_rot=false;};
  virtual void setNonDynamicLinear() {m_is_dynamic=false;};
  virtual void setNonDynamicRot() {m_is_rot=false;};

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
  void resetRotation(){m_quat=Quaternion(1.0,Vec3(0.0,0.0,0.0));}

  friend ostream& operator<<(ostream&, const CRotParticle&);
  void print(){cout << *this << endl << flush;};

  // -- checkpointing --
  virtual void saveSnapShotData(std::ostream& oStream);
  virtual void saveCheckPointData(std::ostream& oStream);
  virtual void loadCheckPointData(std::istream& iStream);

  CRotParticle::exchangeType getExchangeValues();
  void setExchangeValues(const CRotParticle::exchangeType &e);

  template <typename TmplVisitor>
  void visit(TmplVisitor &visitor)
  {
    visitor.visitRotParticle(*this);
  }

  // stress 
  inline double sigma_xx_2D() const {return m_sigma(0,0)/(M_PI*m_rad*m_rad);};
  inline double sigma_xy_2D() const {return m_sigma(0,1)/(M_PI*m_rad*m_rad);};
  inline double sigma_yy_2D() const {return m_sigma(1,1)/(M_PI*m_rad*m_rad);};
//  inline double sigma_d() const;
  static void get_type() {cout <<" CRotParticle" ;};
  friend class TML_PackedMessageInterface;
};

#endif //__ROTPARTICLE_H
