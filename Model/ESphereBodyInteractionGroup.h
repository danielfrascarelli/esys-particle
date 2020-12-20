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

#ifndef __ESPHEREBODYINTERACTIONGROUP_H
#define __ESPHEREBODYINTERACTIONGROUP_H

//--- project includes ---
#include "Model/ESphereBodyInteraction.h"
#include "Model/SphereBodyIG.h"
#include "Model/ElasticInteractionGroup.h"
#include "Model/IGParam.h"

template <class T> class ParallelParticleArray;

//--- STL includes ---
#include <map>

using std::map;

/*!
  \brief Interaction group parameters for CESphereBodyInteractionGroups
*/
class CESphereBodyIGP : public CElasticIGP
{
protected:
  std::string m_spherename;
public:

  CESphereBodyIGP(const std::string&,const std::string&,double);
  virtual void packInto(CVarMPIBuffer*) const;
  std::string getSphereBodyName() const {return m_spherename;};
  friend ostream& operator<<(ostream&,const CESphereBodyIGP&);
};

CESphereBodyIGP* extractESphereBodyIGP(AMPIBuffer*);


// --- Forward decl ---
template <class T> class CESphereBodyInteractionGroup;
template <class T> ostream& operator<< (ostream &, const CESphereBodyInteractionGroup<T> &);

/*!
  \brief Class for a group of unbonded,elastic interactions between particles and a sphere body
*/
template<class T>
class CESphereBodyInteractionGroup : public ASphereBodyInteractionGroup<T>
{
 protected:
  vector<CElasticSphereBodyInteraction<T> > m_interactions;
  double m_k; //!< Elastic modulus
  double m_k_global; //!< total sphere body stiffness
	double k_local;

 public:
  CESphereBodyInteractionGroup(TML_Comm*);
  CESphereBodyInteractionGroup(TML_Comm*,CSphereBody*,const CESphereBodyIGP*);
  virtual ~CESphereBodyInteractionGroup(){}

  /**
   * Null op, time step size not required.
   */
  virtual void setTimeStepSize(double dt)
  {
  }
  
  virtual void calcForces();
  virtual void applyForce(const Vec3&);
  virtual void Update(ParallelParticleArray<T>*);

  friend ostream& operator<< <>(ostream &, const CESphereBodyInteractionGroup &);
};

#include  "ESphereBodyInteractionGroup.hpp"

#endif //__ESPHEREBODYINTERACTIONGROUP_H
