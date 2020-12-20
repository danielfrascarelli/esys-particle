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

#include "Geometry/LayeredBlock.h"

CLayeredBlock2D::CLayeredBlock2D(double xmin,double xmax,double ymin,double ymax,double rmin,double rmax):CRandomBlock2D(xmin,xmax,ymin,ymax,rmin,rmax,1.05)
{}

CLayeredBlock2D::~CLayeredBlock2D()
{}

void CLayeredBlock2D::addLayerBoundary(double d)
{
  LayerBoundaries.insert(d);;
}

void CLayeredBlock2D::generate(int tries,unsigned int seed)
{
  // generate particles
  CRandomBlock2D::generate(tries,seed);
  //-- set tags according to layer --
  int nlayer=0;
  for(set<double>::iterator it1=LayerBoundaries.begin();
      it1!=LayerBoundaries.end();
      it1++){
    nlayer++;
    cout << "layer "<< nlayer << " bdry: " << *it1 << endl;
   for(vector<SimpleParticle>::iterator iter=m_bpart.begin();
	iter!=m_bpart.end();
	iter++){
     if(iter->getPos().Y()>*it1){
	iter->setTag(nlayer);
      }
    }
  }
}
