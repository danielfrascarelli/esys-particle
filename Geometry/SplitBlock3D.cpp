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

#include "Geometry/SplitBlock3D.h"
#include <cstdlib>

/*!
  Constructor of CRandomBlock

  \param xmin minimum in x-direction
  \param xmax maximum in x-direction
  \param ymin minimum in y-direction
  \param ymax maximum in y-direction
  \param zmin minimum in z-direction
  \param zmax maximum in z-direction
  \param rmin minimum particle radius
  \param rmax maximum particle radius
  \param ysplit the position of the split plane 
  \param dir the direction of the split plane (2=y,3=z)
  \param circ_x circular or open boundary conditions in x-direction
*/
CSplitBlock3D::CSplitBlock3D(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double rmin,double rmax,double ysplit,int dir, bool circ_x,bool rough):CRandomBlock3D(xmin,xmax,ymin,ymax,zmin,zmax,rmin,rmax,1.05,circ_x)
{
  m_ysplit=ysplit;
  m_dir=dir;
  if(!rough){
    switch(m_dir){
    case 2 : // split plane y
      { 
	Borders.push_back(Plane3D(Vec3(0.0,1.0,0.0),Vec3(0.0,m_ysplit,0.0))); 
	cout << "split plane y" << endl;
      }
      break;
    case 3 : // split plane z
      { 
	Borders.push_back(Plane3D(Vec3(0.0,0.0,1.0),Vec3(0.0,0.0,m_ysplit))); 
	cout << "split plane z" << endl;
      }
      break;
    default:
      {
	cerr << "invalid direction " << dir << " in CSplitBlock3D" << endl;
      }
    }
  }
}

CSplitBlock3D::~CSplitBlock3D()
{}

/*!
  Fill the space in the block

  \param tries number of times the insertion of a particle is tried 
  \param seed seed for the random number generator
*/
void CSplitBlock3D::generate(int tries,unsigned int seed)
{
    srand(seed);
  // get limits 
  int imin=int(floor(m_xmin/(m_rmax*2.0)));
  int jmin=int(floor(m_ymin/(m_rmax*sqrt(3.0))));
  int kmin=int(floor(m_zmin/(m_rmax*2.0*sqrt(2.0/3.0))));
  int imax=int(ceil(m_xmax/(m_rmax*2.0)));
  int jmax=int(ceil(m_ymax/(m_rmax*sqrt(3.0))));
  int kmax=int(ceil(m_zmax/(m_rmax*2.0*sqrt(2.0/3.0))));
  // do the seeding 
  for(int i=imin;i<=imax;i++){
    for(int j=jmin;j<=jmax;j++){
      for(int k=kmin;k<kmax;k++){
        // calc random radius
        double r=m_random(m_rmin,m_rmax);
        // get position
        double px=(double(i)+0.5*double(j%2)+0.5*double(k%2))*m_rmax*2.0;
        double py=(double(j)+double(k%2)/3.0)*sqrt(3.0)*m_rmax;
        double pz=(double(k)*2.0*sqrt(2.0/3.0))*m_rmax;
        SimpleParticle Po=SimpleParticle(Vec3(px,py,pz),r,getNParts());
        bool fit=checkAFit(Po);
        if(fit){
          insertParticle(Po);	
        }
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
    double p1=0.0;
    double p2=0.0;
    switch(m_dir){
      case 2:
      {
        p1=m_bpart[iter->first()].getPos().Y();
        p2=m_bpart[iter->second()].getPos().Y();
        break;
      }
      case 3:
      {
        p1=m_bpart[iter->first()].getPos().Z();
        p2=m_bpart[iter->second()].getPos().Z();
        break;
      }
      default: break;
    }
    double d1=m_ysplit-p1;
    double d2=m_ysplit-p2;
    if(d1*d2<0){ // on different sides
      set<BasicInteraction,BILess>::iterator h=iter;
      iter++;
      m_iset.erase(h);
    }
  }
}

/*!
  Tag particles along the split plane

  \param tag1 the tag for particles "above" the split (y>y_split)
  \param tag2 the tag for particles "below" the split (y<y_split)
  \param d maximum distance from the split plane at which a particle gets tagged
*/
void CSplitBlock3D::tagSplit(int tag1,int tag2,double d)
{
  for(vector<SimpleParticle>::iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
    double p=0.0;
    switch(m_dir){
    case 2: p=iter->getPos().Y(); break;
    case 3: p=iter->getPos().Z(); break;
    }
    double r=iter->getRad();
    double dist=p-m_ysplit;
    if(fabs(dist)<r+d){ // close to split
      if(dist>0) { // above
	iter->setTag(tag1);
      } else { // below
	iter->setTag(tag2);
      }
    }
  }
}
