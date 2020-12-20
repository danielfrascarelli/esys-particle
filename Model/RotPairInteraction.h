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

#ifndef __ROTPAIRINTERCTION_H
#define __ROTPAIRINTERCTION_H

// -- project includes --
#include "Model/RotParticle.h"
#include "Model/Interaction.h"

/*!
  \class ARotPairInteraction
  \brief Abstract base class interactions between 2 rotational particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
class ARotPairInteraction : public AInteraction
{
 protected:
  CRotParticle *m_p1,*m_p2;

 public:
  // functions 
  ARotPairInteraction();
  ARotPairInteraction(CRotParticle*,CRotParticle*);
  virtual ~ARotPairInteraction();

  inline const CParticle* first()const {return m_p1;}; 
  inline const CParticle* second()const {return m_p2;};
  inline CRotParticle* first() {return m_p1;};
  inline CRotParticle* second() {return m_p2;};
  inline pair<int,int> getPairID() const {return make_pair(m_p1->getID(),m_p2->getID());}
  virtual Vec3 getPos() const = 0;
  virtual void calcForces()=0;
  void checkIDs();
  virtual bool hasTag(int,int) const;
  virtual Vec3 getPosFirst() const {return m_p1->getPos();}
  virtual Vec3 getPosSecond() const {return m_p2->getPos();}
  void setPP(CRotParticle*,CRotParticle*);
  void setPP(const vector<CRotParticle*>);
  
  esys::lsm::quintuple<Vec3,double,Vec3,double,Vec3> getRaw2Data() const
  {
    return 
      esys::lsm::quintuple<Vec3,double,Vec3,double,Vec3>(
        m_p1->getPos(),
        m_p1->getRad(),
        m_p2->getPos(),
        m_p2->getRad(),
        getPos()
      );
  }

  virtual void calcHeatTrans() {}
  virtual void calcHeatFrict() {}

  // dummy implementations for save/load of restart parameters
  virtual void saveRestartData(std::ostream &oStream){};
  virtual void loadRestartData(std::istream &iStream){};
};

#endif // __ROTPAIRINTERCTION_H
