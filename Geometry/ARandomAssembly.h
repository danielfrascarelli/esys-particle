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

#ifndef __ARANDOMASSEMBLY2D_H
#define __ARANDOMASSEMBLY2D_H

//-- project includes --
#include "Geometry/SimpleParticle.h"
#include "Geometry/BasicInteraction.h"
#include "Geometry/SimpleNTable.h"
#include "Geometry/Sphere2d.h"
#include "Geometry/Line.h"

//-- STL includes --
#include <vector>
#include <set>

using std::set;
using std::vector; 

//--- IO includes ---
#include <iostream>

/*!
  \class ARandomAssembly
  \brief Abstract base class for random assemblies, to be used for initialization of random lattices. 
    
  \author Steffen Abe
  $Revision$
  $Date$
*/
class ARandomAssembly
{
 protected:
  ASimpleNTable *m_snt;
  static double m_small_value;
  set<BasicInteraction,BILess> m_iset;
  vector<SimpleParticle> m_bpart;
 
  double m_random(double,double);
  vector<SimpleParticle> getNeighborList(const SimpleParticle&);
  vector<SimpleParticle> get3ClosestNeighbors(const SimpleParticle&, const vector<SimpleParticle>&);
  vector<SimpleParticle> getClosestNeighbors(const SimpleParticle&, int);
  SimpleParticle getClosestParticle(const SimpleParticle&, const vector<SimpleParticle>&);
 
 public:
  virtual ~ARandomAssembly()
  {
  }

  virtual void generate(int,unsigned int)=0;
  virtual void insertParticle(const SimpleParticle)=0;
  virtual void tagParticleClosestTo(const Vec3&,int)=0; 
  virtual void tagEdgeY(int,int,double)=0;
  virtual void tagEdgeZ(int,int,double)=0;
  virtual void tagSplit(int,int,double){std::cout <<"ARA::tagSplit" << std::endl;};

  virtual void writeToGeoFile(const string&)=0;
  virtual void writeToVtkFile(const string&); //{std::cerr << "writeToVtkFile not implemented" << std::endl;}; // empty default implementation
  virtual double calcPorosity()=0;
  virtual vector<pair<double,double> > getSizeDistribution(int)=0;
};

#endif // __ARANDOMASSEMBLY2D_H
