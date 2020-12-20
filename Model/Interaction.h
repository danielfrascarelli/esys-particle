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

#ifndef __INTERACTION_H
#define __INTERACTION_H

// -- project includes --
#include "Model/Particle.h"
#include "Foundation/quintuple.h"

// -- STL includes --
#include <utility>
#include <vector>
using std::vector;
using std::pair;
using std::make_pair;


/*!
  \class AInteraction
  \brief Abstract base class for interactions
  \author Steffen Abe

  $Revision$
  $Date$
*/
class AInteraction
{
 protected:
  vector<int> m_id; //!< id's of the particles involved
  int m_iid;  //!< interaction id
  bool m_init;

 public:
  AInteraction();
  virtual ~AInteraction();

  bool initialized() const;
  virtual void calcForces()=0;
  virtual void calcHeatTrans() {};
  virtual void calcHeatFrict() {};
  vector<int> getAllID() const;
  int getID(){return m_iid;};
  virtual bool hasTag(int,int) const=0;
  virtual Vec3 getPosFirst() const=0;
  virtual Vec3 getPosSecond() const{return Vec3(0.0,0.0,0.0);};
  inline double Count() const {return 1.0;};
};


/*!
  \class APairInteraction
  \brief Abstract base class for 2-particle interactions

  \author Steffen Abe
  $Revision$
  $Date$
*/
class APairInteraction : public AInteraction
{
 protected:
  CParticle *m_p1,*m_p2;

 public:
  // functions 
  APairInteraction();
  APairInteraction(CParticle*,CParticle*);
  virtual ~APairInteraction();
  
  inline const CParticle* first()const {return m_p1;}
  inline const CParticle* second()const {return m_p2;}
  inline CParticle* first() {return m_p1;}
  inline CParticle* second() {return m_p2;}
  
  inline pair<int,int> getPairID() const {return make_pair(m_p1->getID(),m_p2->getID());};
  virtual Vec3 getPos() const = 0;
  virtual void calcForces() = 0;
  void setPP(CParticle*,CParticle*);
  void checkIDs();
  virtual bool hasTag(int,int) const;
  virtual Vec3 getPosFirst() const {return m_p1->getPos();};
  virtual Vec3 getPosSecond() const{return m_p2->getPos();};

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

  template <class TmplParticle> void setPP(const vector<TmplParticle *> &pp)
  {
    m_p1=pp[0];
    m_p2=pp[1];
    m_id.clear();
    m_id.push_back(m_p1->getID());
    m_id.push_back(m_p2->getID());
  }

  // dummy implementations for save/load of restart parameters
  virtual void saveRestartData(std::ostream&){};
  virtual void loadRestartData(std::istream&){};
};
#endif
