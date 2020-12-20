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

#ifndef __VISCWALLIG_H
#define __VISCWALLIG_H

//--- project includes ---
#include "Model/ViscWallInteraction.h"
#include "Model/EWallInteraction.h"
#include "Model/WallIG.h"
#include "Model/EWallInteractionGroup.h"

//--- STL includes ---
#include <map>

using std::map;

template <class T> class ParallelParticleArray;

/*!
  \brief Interaction group parameters for CBWallInteractionGroups
*/
class CVWallIGP : public CEWallIGP
{
 protected:
  int m_tag;
  double m_nu;

 public:
  CVWallIGP(const string&,const string&,double,double,int);
  virtual void  packInto(CVarMPIBuffer*) const;
  void setTag(int tag){m_tag=tag;};
  int getTag()const{return m_tag;};
  void setNu(double nu){m_nu=nu;};
  double getNu()const{return m_nu;};

  friend ostream& operator<<(ostream&,const CVWallIGP&);
};

CVWallIGP* extractVWallIGP(AMPIBuffer*);

// --- Forward decl ---
template <class T> class CViscWallIG;
template <class T> ostream& operator<<(ostream &, const CViscWallIG<T> &);

/*!
  \brief Class for a group of viscous and elastic interactions between particles and a wall
*/
template<class T>
class CViscWallIG : public AWallInteractionGroup<T>
{
protected:
  vector<CViscWallInteraction<T> >    m_visc_interactions;    //!< visc interactions for tagged particles
  vector<CElasticWallInteraction<T> > m_elastic_interactions; //!< elastic interactions all particles
  double m_k; //!< spring constant
  double m_nu;
  int m_tag;

public:
  CViscWallIG(TML_Comm*);
  CViscWallIG(TML_Comm*,CWall*,const CVWallIGP*);
  virtual ~CViscWallIG(){}

  /**
   * Null op, don't require time step size.
   */
  virtual void setTimeStepSize(double dt)
  {
  }

  virtual void calcForces();
  virtual void applyForce(const Vec3&);
  virtual void setVelocity(const Vec3&);
  virtual void Update(ParallelParticleArray<T>*);

  friend ostream& operator<< <>(ostream &, const CViscWallIG &);
};

#include  "Model/ViscWallIG.hpp"

#endif // __VISCWALLIG_H
