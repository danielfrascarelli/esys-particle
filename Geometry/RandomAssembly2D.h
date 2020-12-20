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

#ifndef __RANDOMASSEMBLY2D_H
#define __RANDOMASSEMBLY2D_H

//-- project includes --
#include "Geometry/ARandomAssembly.h"
#include "Geometry/SimpleParticle.h"
#include "Geometry/BasicInteraction.h"
#include "Geometry/Sphere2d.h"
#include "Geometry/Line.h"

//-- STL includes --
#include <vector>

using std::vector; 

/*!
  \class ARandomAssembly2D
  \brief Abstract base class for random assemblies, to be used for initialization of random lattices. 
    
  \author Steffen Abe
  $Revision$
  $Date$
*/
class ARandomAssembly2D : public ARandomAssembly
{
 protected:
  vector<Line> Borders;
  double m_rmin,m_rmax;               //!< min/max particle radius
  double m_xmin,m_xmax,m_ymin,m_ymax; //!< x,y borders of the lattice
  bool m_circ_x;


  virtual Vec3 getAPoint()=0;

  bool isInSpace(const Vec3&);
  bool findAFit(SimpleParticle&, const vector<SimpleParticle>&, const Line&);
  bool findAFit(SimpleParticle&, const vector<SimpleParticle>&);
  virtual bool checkAFit(const SimpleParticle&);
  virtual Line *getClosestPlane(const SimpleParticle&);
  void fillSpace(int);
  virtual int getNParts() const=0;

 public:
  virtual void generate(int,unsigned int)=0;
  virtual void insertParticle(const SimpleParticle)=0;
  virtual void tagParticleClosestTo(const Vec3&,int)=0; 
  virtual void tagEdgeY(int,int,double)=0; 
};

#endif // __RANDOMASSEMBLY2D_H
