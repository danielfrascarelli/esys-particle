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

#include "Geometry/RoughPaddedBlock3d.h"
#include "Geometry/RandomAssembly3D.h"

/*!
  Get closest plane

  \param Po the particle
*/
Plane3D CRoughPaddedBlock3D::getClosestPlane(const SimpleParticle& Po)
{
  Plane3D PL=ARandomAssembly3D::getClosestPlane(Po);
  double spl=PL.sep(Po.getPos());

  RectPatch RP=getClosestPatch(Po,spl);
  double srp=RP.sep(Po.getPos());
  if (srp<spl){
    PL=RP.getPlane(Po.getPos());
  }

  return PL;
}


/*!
  Get closest fault patch

  \param Po the particle
*/
RectPatch CRoughPaddedBlock3D::getClosestPatch(const SimpleParticle& Po, double sep0)
{
  RectPatch RP=(*m_fault.begin());
  double oldsep=sep0;

  for(vector<RectPatch>::iterator iter=m_fault.begin();
      iter!=m_fault.end();
      iter++){
    double sep=iter->sep(Po.getPos());
    if((sep!=-1.0) && ((sep<oldsep) || (oldsep==-1.0))) {
      RP=*iter;
      oldsep=sep;
    }
  }
  return RP;
}

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
  \param ysplit
  \param pad
  \param circ_x
*/
CRoughPaddedBlock3D::CRoughPaddedBlock3D(double xmin, double xmax,
					 double ymin ,double ymax,
					 double zmin, double zmax,
					 double rmin, double rmax,
					 double ysplit, double pad, bool circ_x):
  CPaddedBlock3D(xmin,xmax,ymin,ymax,zmin,zmax,rmin,rmax,ysplit,pad,2,circ_x)
{
  std::cout << "CRoughPaddedBlock3D" << std::endl;
  // ugly trick: remove last plane (splitplane)
  Borders.pop_back();
}

/*!
  setup fault roughness 

  \param nx x-resolution
  \param nz z-resolution
  \param depth amount of roughness in the rough patches
  \param prob probability of a patch to be rough
*/
void CRoughPaddedBlock3D::setRoughness(int nx,int nz,double depth,double prob)
{
  double dx=(m_xmax-m_xmin)/double(nx);
  double dz=(m_zmax-m_zmin)/double(nz);

  std::cout << "dx,dz: " << dx << ", " << dz << std::endl;

  for(int i=0;i<nx;i++){
    for(int j=0;j<nz;j++){
      double xmin=double(i)*dx;
      double xmax=double(i+1)*dx;
      double zmin=double(j)*dz;
      double zmax=double(j+1)*dz;
      double dy=(m_random(0.0,1.0)<prob) ? depth : 0.0;
      m_fault.push_back(RectPatch(xmin,xmax,zmin,zmax,m_ysplit,dy));	
      std::cout << dy << " "; 
    }    
    std::cout << std::endl;
  }
}

/*!
  check if Po is within the Space and is not crossing any boundary or 
  overlapping with other particles.

  \param Po the particle
*/
bool CRoughPaddedBlock3D::checkAFit(const SimpleParticle& P) 
{
  bool fit=ARandomAssembly3D::checkAFit(P);

  // check against fault
  vector<RectPatch>::iterator iter=m_fault.begin();
  while((iter!=m_fault.end()) && (fit)){
    double s=iter->dist(P.getPos());
     fit=(s>P.getRad()-1e-4) || (s==-1.0); // remove -1 case once proper dist implemented
     iter++;
   }
  return fit;
}

/*!
  generate the particle packing

  \param tries number of attempts to insert particle before giving up
  \param seed random seed 
*/
//void CRoughPaddedBlock3D::generate(int tries,unsigned int seed)
void CRoughPaddedBlock3D::generate(int tries)
{
  generate_regular_padding();

  // random block
  // limits
  int imin=int(floor(m_xmin/(m_rmax*2.0)));
  int imax=int(ceil(m_xmax/(m_rmax*2.0)));
  int jmin=int(floor(m_zmin/(m_rmax*sqrt(3.0))));
  int jmax=int(ceil(m_zmax/(m_rmax*sqrt(3.0))));
  int kmin=int(floor((m_ymin+m_pad_size-m_rmax)/(m_rmax*2.0*sqrt(2.0/3.0))));
  int kmax=int(ceil((m_ymax-(m_pad_size-m_rmax))/(m_rmax*2.0*sqrt(2.0/3.0))));
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

  // remove interactions crossing any fault segment
  // for(vector<RectPatch>::iterator iter=m_fault.begin();
//       iter!=m_fault.end();
//       iter++){
//     for(set<BasicInteraction,BILess>::iterator iter_in=m_iset.begin();
// 	iter_in!=m_iset.end();
// 	iter_in++){
//       Vec3 P1=m_bpart[iter_in->first()].getPos();
//       Vec3 P2=m_bpart[iter_in->second()].getPos();
//       if(iter->intersect(P1,P2)){ // on different sides
// 	set<BasicInteraction,BILess>::iterator h=iter_in;
// 	iter_in++;
// 	m_iset.erase(h);
// 	iter_in--;
//       }
//     }   
//   }
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

// void CRoughPaddedBlock3D::tagSplit(int,int,double)
// {
//   CSplitBlock(
// }
