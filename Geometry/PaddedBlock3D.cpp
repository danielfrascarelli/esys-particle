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

#include "Geometry/PaddedBlock3D.h"
#include <cstdlib>

/*!
  Constructor of CPaddedBlock3d

  \param xmin minimum in x-direction
  \param xmax maximum in x-direction
  \param ymin minimum in y-direction
  \param ymax maximum in y-direction
  \param zmin minimum in z-direction
  \param zmax maximum in z-direction
  \param rmin minimum particle radius
  \param rmax maximum particle radius
*/
CPaddedBlock3D::CPaddedBlock3D(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double rmin,double rmax,double ysplit,double pad,int dir,bool circ_x):CSplitBlock3D(xmin,xmax,ymin,ymax,zmin,zmax,rmin,rmax,ysplit,dir,circ_x)
{
  m_pad_size=pad;
}

/*!
  Generate a random point within the space of the random part of the block
*/

Vec3 CPaddedBlock3D::getAPoint()
{
  double px=0.0;
  double py=0.0;
  double pz=0.0;
  
  px=m_random(m_xmin+m_rmin,m_xmax-m_rmin);
  switch(m_dir){
  case 2:
    {
      py=m_random(m_ymin+m_pad_size-m_rmin,m_ymax-m_pad_size+m_rmin);
      pz=m_random(m_zmin+m_rmin,m_zmax-m_rmin);
    }
    break;
  case 3:
    {
      py=m_random(m_ymin+m_rmin,m_ymax-m_rmin);
      pz=m_random(m_zmin+m_pad_size-m_rmin,m_zmax-m_pad_size+m_rmin);
    }
    break;
  }
  return Vec3(px,py,pz);
}

/*!
  generate regular padding sections
*/
void CPaddedBlock3D::generate_regular_padding()
{
  // lower padding area
  // limits
  int kmin=0;
  int kmax=0;
  int jmin=0;
  int jmax=0;
  // x-dir
  int imin=int(floor(m_xmin/(m_rmax*2.0)));
  int imax=int(ceil((m_xmax+m_rmax)/(m_rmax*2.0)));
  switch(m_dir){
  case 2:
    {
      jmin=int(floor(m_zmin/(m_rmax*sqrt(3.0))));
      kmin=int(floor(m_ymin/(m_rmax*2.0*sqrt(2.0/3.0))));
      jmax=int(ceil(m_zmax/(m_rmax*sqrt(3.0))));
      kmax=int(ceil((m_ymin+m_pad_size-m_rmax)/(m_rmax*2.0*sqrt(2.0/3.0))));
    }
    break;
  case 3:
    {
      jmin=int(floor(m_zmin/(m_rmax*sqrt(3.0))));
      kmin=int(floor(m_ymin/(m_rmax*2.0*sqrt(2.0/3.0))));
      jmax=int(ceil((m_zmin+m_pad_size-m_rmax)/(m_rmax*sqrt(3.0))));
      kmax=int(ceil(m_ymax/(m_rmax*2.0*sqrt(2.0/3.0))));
    }
    break;
  }
  // particles
  // lower padding area
  for(int i=imin;i<=imax;i++){
    for(int j=jmin;j<jmax;j++){
      for(int k=kmin;k<kmax;k++){
        double px=(double(i)+0.5*double(j%2)+0.5*double(k%2))*m_rmax*2.0;
        double pz=((double(j)+double(k%2)/3.0)*sqrt(3.0)+1.0)*m_rmax;
        double py=(double(k)*2.0*sqrt(2.0/3.0)+1.0)*m_rmax;
        SimpleParticle Po=SimpleParticle(Vec3(px,py,pz),m_rmax,getNParts());
        bool fit=checkAFit(Po);
        if(fit){
          insertParticle(Po);	
        }
      }
    }
  }
  // upper padding area
  switch(m_dir){
  case 2:
    {
      for(int i=imin;i<=imax;i++){
        for(int j=jmin;j<jmax;j++){
          for(int k=kmin;k<kmax;k++){
            double px=(double(i)+0.5*double(j%2)+0.5*double(k%2))*m_rmax*2.0;
            double pz=((double(j)+double(k%2)/3.0)*sqrt(3.0)+1.0)*m_rmax;
            double py=(m_ymax-m_ymin)-(double(k)*2.0*sqrt(2.0/3.0)+1.0)*m_rmax;
            SimpleParticle Po=SimpleParticle(Vec3(px,py,pz),m_rmax,getNParts());
            bool fit=checkAFit(Po);
            if(fit){
              insertParticle(Po);	
            }
          }
        }
      }
    }
    break;
  case 3:
    {
      for(int i=imin;i<=imax;i++){
        for(int j=jmin;j<jmax;j++){
          for(int k=kmin;k<kmax;k++){
            double px=(double(i)+0.5*double(j%2)+0.5*double(k%2))*m_rmax*2.0;
            double pz=(m_zmax-m_zmin)-((double(j)+double(k%2)/3.0)*sqrt(3.0)+1.0)*m_rmax;
            double py=(double(k)*2.0*sqrt(2.0/3.0)+1.0)*m_rmax;
            SimpleParticle Po=SimpleParticle(Vec3(px,py,pz),m_rmax,getNParts());
            bool fit=checkAFit(Po);
            if(fit){
              insertParticle(Po);	
            }
          }
        }
      }
    }
    break;
  }

}

/*!
  Fill the space in the block

  \param tries number of times the insertion of a particle is tried 
  \param seed seed for the random number generator
*/
void CPaddedBlock3D::generate(int tries,unsigned int seed)
{
  srand(seed);
  
  generate_regular_padding();


  // random block
  // limits
  int jmin,jmax,kmin,kmax;

  int imin=int(floor(m_xmin/(m_rmax*2.0)));
  int imax=int(ceil(m_xmax/(m_rmax*2.0)));
  switch(m_dir){
  case 2:
    {
      jmin=int(floor(m_zmin/(m_rmax*sqrt(3.0))));
      jmax=int(ceil(m_zmax/(m_rmax*sqrt(3.0))));
      kmin=int(floor((m_ymin+m_pad_size-m_rmax)/(m_rmax*2.0*sqrt(2.0/3.0))));
      kmax=int(ceil((m_ymax-(m_pad_size-m_rmax))/(m_rmax*2.0*sqrt(2.0/3.0))));
    }
    break;
  case 3:
    {
      jmin=int(floor((m_zmin+m_pad_size-m_rmax)/(m_rmax*sqrt(3.0))));
      jmax=int(ceil((m_zmax-(m_pad_size-m_rmax))/(m_rmax*sqrt(3.0))));
      kmin=int(floor(m_ymin/(m_rmax*2.0*sqrt(2.0/3.0))));
      kmax=int(ceil(m_ymax/(m_rmax*2.0*sqrt(2.0/3.0))));
    } 
    break;
  default:
    {
      // should never get here
      jmin=jmax=kmin=kmax=0;
    }
  }
  // particles
  for(int i=imin;i<=imax;i++){
    for(int j=jmin;j<jmax;j++){
      for(int k=kmin;k<kmax;k++){
        // calc random radius
        double r=m_random(m_rmin,m_rmax);
        // get position
        double px=(double(i)+0.5*double(j%2)+0.5*double(k%2))*m_rmax*2.0;
        double pz=((double(j)+double(k%2)/3.0)*sqrt(3.0))*m_rmax;
        double py=(double(k)*2.0*sqrt(2.0/3.0)+1.0)*m_rmax;
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
      iter--;
    } 
  }
}  
