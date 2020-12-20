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

#ifndef __SOFTBWALLINTERACTIONGROUP_H
#define __SOFTBWALLINTERACTIONGROUP_H

//--- project includes ---
#include "Model/SoftBWallInteraction.h"
#include "Model/WallIG.h"
#include "Model/BWallInteractionGroup.h"

//--- STL includes ---
#include <map>

using std::map;

template <class T> class ParallelParticleArray;

/*!
  \brief Interaction group parameters for CSoftBWallInteractionGroups
*/
class CSoftBWallIGP : public CBWallIGP
{
 protected:
  double m_shearK;
  bool m_scaling;
 public:
  CSoftBWallIGP(const std::string&,const std::string&,double,double,int,int,bool);
  virtual void  packInto(CVarMPIBuffer*) const;
  double getNormalK()const{return m_k;};
  double getShearK()const{return m_shearK;};
  bool getScaling()const{return m_scaling;};

  friend ostream& operator<<(ostream&,const CSoftBWallIGP&);
};

CSoftBWallIGP* extractSoftBWallIGP(AMPIBuffer*);

// --- Forward decl ---
template <class T> class CSoftBWallInteractionGroup;
template <class T> ostream& operator<<(ostream &, const CSoftBWallInteractionGroup<T> &);


/*!
  \brief Class for a group of bonded, elastic interactions with per-direction spring constants
	between particles and a wall
*/
template<class T>
class CSoftBWallInteractionGroup : public AWallInteractionGroup<T>
{
 protected:
  vector<CSoftBondedWallInteraction<T> > m_interactions;
  double m_normalK,m_shearK; //!< spring constants for each direction
  int m_tag;
  int m_mask;
  bool m_scaling;

 public:
  CSoftBWallInteractionGroup(TML_Comm*);
  CSoftBWallInteractionGroup(TML_Comm*,CWall*,const CSoftBWallIGP*);
  virtual ~CSoftBWallInteractionGroup(){}

  /**
   * Null op, don't require time step size.
   */
  virtual void setTimeStepSize(double dt)
  {
  }

  virtual void calcForces();
  virtual void applyForce(const Vec3&);
  virtual void Update(ParallelParticleArray<T>*);

  friend ostream& operator<< <>(ostream &, const CSoftBWallInteractionGroup &);
};

#include  "SoftBWallInteractionGroup.hpp"

#endif //__BSOFTWALLINTERACTIONGROUP_H
