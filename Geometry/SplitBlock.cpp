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

#include "SplitBlock.h"
#include <cstdlib>

CSplitBlock2D::CSplitBlock2D(double xmin,double xmax,double ymin,double ymax,double rmin,double rmax,double ysplit,bool circ_x):CRandomBlock2D(xmin,xmax,ymin,ymax,rmin,rmax,circ_x)
{
  cout << "CSplitBlock2D" << endl;
  m_ysplit=ysplit;
  Borders.push_back(Line(Vec3(0.0,1.0,0.0),Vec3(0.0,m_ysplit,0.0))); // split line
}

CSplitBlock2D::~CSplitBlock2D()
{}

void CSplitBlock2D::generate(int tries,unsigned int seed)
{
  srand(seed);
  // get limits 
  int imin=int(floor(m_xmin/(m_rmax*2.0)));
  int jmin=int(floor(m_ymin/(m_rmax*sqrt(3.0))));
  int imax=int(ceil(m_xmax/(m_rmax*2.0)));
  int jmax=int(ceil(m_ymax/(m_rmax*sqrt(3.0))));
  // do the seeding 
  for(int i=imin;i<=imax;i++){
    for(int j=jmin;j<=jmax;j++){
      // calc random radius
      double r=m_random(m_rmin,m_rmax);
      // get position
      double px=(double(i)+0.5*double(j%2))*m_rmax*2.0;
      double py=double(j)*sqrt(3.0)*m_rmax;
      SimpleParticle Po=SimpleParticle(Vec3(px,py,0),r,getNParts());
      bool fit=checkAFit(Po);
      if(fit){
	insertParticle(Po);	
      }
    }
  }
  // fill space
  fillSpace(tries);

  // get set of interactions
  m_snt->getInteractions(m_iset,1.05);

  // remove interactions crossing the split
  for(set<BasicInteraction,BILess>::iterator iter=m_iset.begin();
      iter!=m_iset.end();
      iter++){
    double y1=m_bpart[iter->first()].getPos().Y();
    double y2=m_bpart[iter->second()].getPos().Y();
    double dy1=m_ysplit-y1;
    double dy2=m_ysplit-y2;
    if(dy1*dy2<0){ // on different sides
      set<BasicInteraction,BILess>::iterator h=iter;
      iter++;
      m_iset.erase(h);
    }
  }
}

/*!
  Tag particles along the split line

  \param tag1 the tag for particles "above" the split (y>y_split)
  \param tag2 the tag for particles "below" the split (y<y_split)
  \param d maximum distance from the split line at which a particle gets tagged
*/
void CSplitBlock2D::tagSplit(int tag1,int tag2,double d)
{
  for(vector<SimpleParticle>::iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
    double py=iter->getPos().Y();
    double r=iter->getRad();
    double dist=py-m_ysplit;
    if(fabs(dist)<r+d){ // close to split
      if(dist>0) { // above
	iter->setTag(tag1);
      } else { // below
	iter->setTag(tag2);
      }
    }
  }
}
