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

#ifndef __RANDOMASSEMBLY3D_H
#define __RANDOMASSEMBLY3D_H

//-- project includes --
#include "Geometry/ARandomAssembly.h"
#include "Geometry/SimpleParticle.h"
#include "Geometry/BasicInteraction.h"
#include "Geometry/Plane3D.h"

//-- STL includes --
#include <vector>
#include <set>

using std::set;
using std::vector; 

/*!
  \class ARandomAssembly3D
  \brief Abstract base class for random assemblies, to be used for initialization of random lattices. 
    
  \author Steffen Abe
  $Revision$
  $Date$
*/
class ARandomAssembly3D : public ARandomAssembly
{
 protected:
  vector<Plane3D> Borders;
  double m_rmin,m_rmax;                             //!< min/max particle radius
  double m_xmin,m_xmax,m_ymin,m_ymax,m_zmin,m_zmax; //!< x,y,z borders of the lattice
  bool m_circ_x;

  virtual Vec3 getAPoint()=0;

  bool findAFit(SimpleParticle&, const vector<SimpleParticle>&);
  bool findAFit(SimpleParticle&, const vector<SimpleParticle>&, const Plane3D&);
  virtual bool checkAFit(const SimpleParticle&);
  virtual Plane3D getClosestPlane(const SimpleParticle&);
  void fillSpace(int);
  virtual int getNParts() const=0;

 public:
  virtual void generate(int,unsigned int)=0;
  virtual void insertParticle(const SimpleParticle)=0;
  virtual void tagParticleClosestTo(const Vec3&,int)=0; 
  virtual void tagEdgeY(int,int,double)=0; 
  virtual void tagEdgeZ(int,int,double)=0; 
};

#endif // __RANDOMASSEMBLY3D_H
