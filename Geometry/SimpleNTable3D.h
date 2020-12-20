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

#ifndef __SIMPLENTABLE3D_H
#define __SIMPLENTABLE3D_H

//-- project includes --
#include "Geometry/SimpleParticle.h"
#include "Geometry/BasicInteraction.h"
#include "Geometry/SimpleNTable.h"

typedef std::set<BasicInteraction,BILess> InteractionSet;
typedef std::vector<SimpleParticle>       ParticleVector;

/*!
  \class CSimple2DNTable
  \brief 2D implementation of simple, serial neighbor table

  \author Steffen Abe
  $Revision$
  $Date$
*/
class CSimple3DNTable : public ASimpleNTable
{
private:
  Vec3 m_xshift,m_yshift,m_zshift;
  int m_xsize,m_ysize,m_zsize;
  bool m_xcirc,m_ycirc,m_zcirc;

protected:
  virtual int index(const Vec3&) const;
  virtual vector<int> allidx(const Vec3&) const;
  virtual void insertParticleCircular(SimpleParticle);

public:
  CSimple3DNTable(const Vec3&,const Vec3&,double,bool xcirc=false,bool ycirc=false,bool zcirc=false);
  virtual void getInteractions(set<BasicInteraction,BILess>&,double);
  void print();
  
  template <class TmplInteractionValidator>
  InteractionSet getInteractions(const TmplInteractionValidator &validator) const
  {
    InteractionSet iset;
    for(int i=0;i<m_xsize;i++){
      for(int j=0;j<m_ysize;j++){
        for(int k=0;k<m_zsize;k++){
          int idx=i+m_xsize*j+k*m_xsize*m_zsize;
          if(m_data[idx].size() >= 2){
            for(ParticleVector::const_iterator iter = m_data[idx].begin();
                iter != m_data[idx].end()-1;
                iter++)
            {
              for (
                ParticleVector::const_iterator iter2 = iter+1;
                iter2 != m_data[idx].end();
                iter2++)
              {
                if (validator.isValid(*iter, *iter2))
                {
                  iset.insert(BasicInteraction(iter->getID(),iter2->getID()));
                }
              }
            }
          }
        }
      }
    }
    return iset;
  }
};

#endif //__SIMPLENTABLE3D_H
