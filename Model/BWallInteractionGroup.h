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


#ifndef __BWALLINTERACTIONGROUP_H
#define __BWALLINTERACTIONGROUP_H

//--- project includes ---
#include "Model/BWallInteraction.h"
#include "Model/EWallInteraction.h"
#include "Model/WallIG.h"
#include "Model/EWallInteractionGroup.h"

//--- STL includes ---
#include <map>

using std::map;

template <class T> class ParallelParticleArray;

/*!
  \class CBWallIGP
  \brief Interaction group parameters for CBWallInteractionGroups

  \author Steffen Abe
  $Revision$
  $Date$
*/
class CBWallIGP : public CEWallIGP
{
 protected:
  int m_tag;
  int m_mask;

 public:
  CBWallIGP(const std::string&,const std::string&,double,int,int);
  virtual void  packInto(CVarMPIBuffer*) const;
  int getTag()const{return m_tag;};
  int getMask()const{return m_mask;};

  friend ostream& operator<<(ostream&,const CBWallIGP&);
};

CBWallIGP* extractBWallIGP(AMPIBuffer*);

// --- Forward decl ---
template <class T> class CBWallInteractionGroup;
template <class T> ostream& operator<<(ostream &, const CBWallInteractionGroup<T> &);

/*!
  \class CBWallInteractionGroup
  \brief Class for a group of bonded,elastic interactions between particles and a wall

  \author Steffen Abe
  $Revision$.
  $Date$
*/
template<class T>
class CBWallInteractionGroup : public AWallInteractionGroup<T>
{
 protected:
  vector<CBondedWallInteraction<T> > m_bonded_interactions; //!< bonded interactions for tagged particles
  vector<CElasticWallInteraction<T> > m_elastic_interactions; //!< elastic interactions for the rest
  double m_k; //!< spring constant
  int m_tag;
  int m_mask;
  
 public:
  CBWallInteractionGroup(TML_Comm*);
  CBWallInteractionGroup(TML_Comm*,CWall*,const CBWallIGP*);
  virtual ~CBWallInteractionGroup(){}

  virtual void calcForces();
  virtual void applyForce(const Vec3&);
  virtual void Update(ParallelParticleArray<T>*);

  friend ostream& operator<< <> (ostream &, const CBWallInteractionGroup &);
};

#include  "BWallInteractionGroup.hpp"

#endif //__BWALLINTERACTIONGROUP_H
