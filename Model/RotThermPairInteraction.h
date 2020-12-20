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

#ifndef __ROTTHERMPAIRINTERCTION_H
#define __ROTTHERMPAIRINTERCTION_H

// -- project includes --
#include "Model/RotThermParticle.h"
#include "Model/Interaction.h"

/*!
  Abstract base class interactions between 2 rotational particles
*/
class ARotThermPairInteraction : public AInteraction
{
 protected:
  CRotThermParticle *m_p1,*m_p2;

 public:
  // functions 
  ARotThermPairInteraction();
  ARotThermPairInteraction(CRotThermParticle*,CRotThermParticle*);
  virtual ~ARotThermPairInteraction();

  inline const CRotThermParticle* first()const {return m_p1;}; 
  inline const CRotThermParticle* second()const {return m_p2;};
  inline CRotThermParticle* first() {return m_p1;};
  inline CRotThermParticle* second() {return m_p2;};
  inline pair<int,int> getPairID() const {return make_pair(m_p1->getID(),m_p2->getID());}
  virtual Vec3 getPos() const = 0;
  virtual void calcForces()=0 ;
  virtual void calcHeatFrict(){} ;
  virtual void calcHeatTrans(){} ;
  void checkIDs();
  virtual bool hasTag(int,int) const;
  virtual Vec3 getPosFirst() const {return m_p1->getPos();}
  virtual Vec3 getPosSecond() const {return m_p2->getPos();}
  void setPP(CRotThermParticle*,CRotThermParticle*);
  void setPP(const vector<CRotThermParticle*>);
  
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
  
};

#endif // __ROTTHERMPAIRINTERCTION_H
