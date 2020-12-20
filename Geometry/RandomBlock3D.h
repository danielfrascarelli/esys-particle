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

#ifndef __RANDOMBLOCK3D_H
#define __RANDOMBLOCK3D_H

//-- project includes --
#include "Geometry/SimpleParticle.h"
#include "Geometry/SimpleNTable.h"
#include "Geometry/Sphere2d.h"
#include "Geometry/Line.h"
#include "Geometry/RandomAssembly3D.h"

//-- STL includes --
#include <vector>
#include <string>
#include <utility>

using std::vector; 
using std::string;
using std::pair;

/*!
  \class CRandomBlock3D
  \brief Class for the generation of a 3D random lattice in a rectangular area. 

  \author Steffen Abe
  $Revision$
  $Date$
*/  
class CRandomBlock3D : public ARandomAssembly3D
{
 private:
  bool m_is_bonded;

 protected:
  virtual Vec3 getAPoint();
  virtual int getNParts() const{return m_bpart.size();};
  double m_maxConnDist;

 public:
  CRandomBlock3D(double,double,double,double,double,double,double,double,double,bool,bool is_bonded=true);
  virtual ~CRandomBlock3D();

  virtual void generate(int,unsigned int);
  virtual void insertParticle(const SimpleParticle);
  //virtual void tagParticleClosestTo(const Vec3&,int); 
  virtual void tagEdgeY(int,int,double);
  virtual void tagEdgeZ(int,int,double);
 
  virtual void writeToGeoFile(const string&);
  virtual double calcPorosity();
  virtual vector<pair<double,double> > getSizeDistribution(int);
};

#endif //  __RANDOMBLOCK3D_H
 
