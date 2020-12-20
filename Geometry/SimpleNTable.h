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

#ifndef __SIMPLENTABLE_H
#define __SIMPLENTABLE_H

//-- project includes --
#include "Geometry/SimpleParticle.h"
#include "Geometry/BasicInteraction.h"

//-- STL-includes --
#include <vector>
#include <set>

using std::vector;
using std::set;

/*!
  \class ASimpleNTable
  \brief Abstract base class providing the interface for a simple, serial neighbor table. Used in random initialization.

  \author Steffen Abe
  $Revision$
  $Date$
*/
class ASimpleNTable
{
protected:
  vector<SimpleParticle> *m_data;
  Vec3                   m_p0;
  double                 m_dim;
  int                    m_numInsertedParticles;  

  virtual int index(const Vec3&) const=0;
  virtual vector<int> allidx(const Vec3&) const=0;
  virtual void insertParticleCircular(SimpleParticle)=0;

public:
  ASimpleNTable();
  virtual ~ASimpleNTable();

  int getNumInsertedParticles() const;
  const vector<SimpleParticle>* getNeighbors(const Vec3&) const;
  int getClosestParticleID(const Vec3&) const;
  virtual void getInteractions(set<BasicInteraction,BILess>&,double)=0;
  void insertParticle(SimpleParticle);
};

/*!
  \class CSimple2DNTable
  \brief 2D implementation of simple, serial neighbor table

  \author Steffen Abe
  $Revision$
  $Date$
*/

class CSimple2DNTable : public ASimpleNTable
{
private:
  Vec3 m_xshift,m_yshift;
  int m_xsize,m_ysize;
  bool m_xcirc,m_ycirc;
  
protected:
  virtual int index(const Vec3&) const;
  virtual vector<int> allidx(const Vec3&) const;
  virtual void insertParticleCircular(SimpleParticle);

public:
  CSimple2DNTable(const Vec3&,const Vec3&,double,bool xcirc=false,bool ycirc=false);
  virtual void getInteractions(set<BasicInteraction,BILess>&,double);
  void print();
};


#endif //__SIMPLENTABLE_H
