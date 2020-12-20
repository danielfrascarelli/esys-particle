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

#ifndef __EWALLINTERACTIONGROUP_H
#define __EWALLINTERACTIONGROUP_H

//--- project includes ---
#include "Model/EWallInteraction.h"
#include "Model/WallIG.h"
#include "Model/ElasticInteractionGroup.h"
#include "Model/IGParam.h"

template <class T> class ParallelParticleArray;

//--- STL includes ---
#include <map>

using std::map;

/*!
  \brief Interaction group parameters for CEWallInteractionGroups
*/
class CEWallIGP : public CElasticIGP
{
protected:
  std::string m_wallname;
public:

  CEWallIGP(const std::string&,const std::string&,double);
  virtual void packInto(CVarMPIBuffer*) const;
  std::string getWallName() const {return m_wallname;};
  friend ostream& operator<<(ostream&,const CEWallIGP&);
};

CEWallIGP* extractEWallIGP(AMPIBuffer*);


// --- Forward decl ---
template <class T> class CEWallInteractionGroup;
template <class T> ostream& operator<< (ostream &, const CEWallInteractionGroup<T> &);

/*!
  \brief Class for a group of unbonded,elastic interactions between particles and a wall
*/
template<class T>
class CEWallInteractionGroup : public AWallInteractionGroup<T>
{
 protected:
  vector<CElasticWallInteraction<T> > m_interactions;
  double m_k; //!< Elastic modulus
  double m_k_global; //!< total wall stiffness

 public:
  CEWallInteractionGroup(TML_Comm*);
  CEWallInteractionGroup(TML_Comm*,CWall*,const CEWallIGP*);
  virtual ~CEWallInteractionGroup(){}

  /**
   * Null op, time step size not required.
   */
  virtual void setTimeStepSize(double dt)
  {
  }
  
  virtual void calcForces();
  virtual void applyForce(const Vec3&);
  virtual void Update(ParallelParticleArray<T>*);

  friend ostream& operator<< <>(ostream &, const CEWallInteractionGroup &);
};

#include  "EWallInteractionGroup.hpp"

#endif //__EWALLINTERACTIONGROUP_H
