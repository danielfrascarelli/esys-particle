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

#include "Geometry/FaultedBlock2d.h"
#include <cstdlib>

// --- IO includes ---
#include "fstream"

using std::ofstream;

/*!
  Constructor. Set up "empty" block, i.e. without any fault segments -> add those via addSegment

  \param xmin minimum in x-direction
  \param xmax maximum in x-direction
  \param ymin minimum in y-direction
  \param ymax maximum in y-direction
  \param rmin minimum particle radius
  \param rmax maximum particle radius
  \param pad thickness of the padding region at each diving edge
  \param circ_x circular boudary condition in x-direction
*/
FaultedBlock2D::FaultedBlock2D(double xmin,double xmax,double ymin,double ymax,double rmin,double rmax,double pad,bool circ_x) :
  CRandomBlock2D(xmin,xmax,ymin,ymax,rmin,rmax,circ_x)
{
  m_pad_size=pad;
}

FaultedBlock2D::~FaultedBlock2D()
{}

/*!
  Get a random point with in the "random" region, i.e. outside the padding zone
*/
Vec3  FaultedBlock2D::getAPoint()
{
  double px,py;
  
  px=m_random(m_xmin+m_rmin,m_xmax-m_rmin);
  py=m_random(m_ymin+m_pad_size+m_rmin,m_ymax-m_pad_size-m_rmin);

  return Vec3(px,py,0.0);
}

/*!
  Add a fault segment to the block

  \param P1 1st end point
  \param P2 2nd end point
  \param r roughness parameter (0.0...1.0)
*/
void FaultedBlock2D::addSegment(const Vec3& P1,const Vec3& P2,double r)
{
  //std::cout << "addSegment " << P1 << " - " << P2 << "  r= " << r << std::endl;
  // calculate normal for shift
  Vec3 U=(P2-P1).unit();
  Vec3 N;
  N.X()=U.Y();
  N.Y()=-U.X();
  N.Z()=0.0;
  Vec3 P1a=P1+0.5*r*N;
  Vec3 P2a=P2+0.5*r*N;
  Vec3 P1b=P1-0.5*r*N;
  Vec3 P2b=P2-0.5*r*N;
  m_fault.push_back(make_pair(r,LineSegment(P1a,P2a)));
  m_fault.push_back(make_pair(r,LineSegment(P2b,P1b)));
  m_f2.push_back(LineSegment(P1,P2));
}

/*!
  check if Po is within the Space and is not crossing any boundary or 
  overlapping with other particles.

  \param Po the particle
*/
bool FaultedBlock2D::checkAFit(const SimpleParticle& Po)
{
  bool fit=ARandomAssembly2D::checkAFit(Po);
  
  if(fit){ // check intersection with fault segments
    for(vector<pair<double,LineSegment> >::iterator iter=m_fault.begin();
	iter!=m_fault.end();
	iter++){
      double d=iter->second.sep(Po.getPos());
      fit=(Po.getRad()<d+iter->first+10e-4);
      //std::cout << "check pos : " << Po.getPos() << " rad " << Po.getRad() << " vs. line segment " << iter->second.getP1() << " - " <<  iter->second.getP1() << " res: " << fit << std::endl;
    }
  }

  return fit;
}

/*!
  generate the particle packing

  \param tries number of attempts to insert particle before giving up
  \param seed random seed 
*/
void FaultedBlock2D::generate(int tries,unsigned int seed)
{
  int imin,imax,jmin,jmax;

  srand(seed);
  imin=int(floor(m_xmin/(m_rmax*2.0)));
  imax=int(ceil(m_xmax/(m_rmax*2.0)));
  // "lower" padding area
  // y-limits
  jmin=int(floor(m_ymin/(m_rmax*sqrt(3.0))));
  jmax=int(floor((m_ymin+m_pad_size)/(m_rmax*sqrt(3.0))))+1;
  // do the seeding 
  for(int i=imin;i<=imax;i++){
    for(int j=jmin;j<=jmax;j++){
      double px=(double(i)+0.5*double(j%2))*m_rmax*2.0;
      double py=double(j)*sqrt(3.0)*m_rmax;
      SimpleParticle Po=SimpleParticle(Vec3(px,py,0),m_rmax,getNParts());
      bool fit=checkAFit(Po);
      if(fit){
         insertParticle(Po);	
      }
    }
  }
  // "upper" padding area
  jmin=int(ceil((m_ymax-m_pad_size)/(m_rmax*sqrt(3.0))))-1;
  jmax=int(ceil(m_ymax/(m_rmax*sqrt(3.0))));
  // do the seeding 
  for(int i=imin;i<=imax;i++){
    for(int j=jmin;j<=jmax;j++){
      double px=(double(i)+0.5*double(j%2))*m_rmax*2.0;
      double py=double(j)*sqrt(3.0)*m_rmax;
      SimpleParticle Po=SimpleParticle(Vec3(px,py,0),m_rmax,getNParts());
      bool fit=checkAFit(Po);
      if(fit){
         insertParticle(Po);	
      }
    }
  }
  // "random" area
  // get limits 
  jmin=int(floor((m_ymin+m_pad_size)/(m_rmax*sqrt(3.0))))+2;
  jmax=int(ceil((m_ymax-m_pad_size)/(m_rmax*sqrt(3.0))))-2;
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

  // remove interactions crossing any fault segment
  for(vector<LineSegment>::iterator iter=m_f2.begin();
      iter!=m_f2.end();
      iter++){
    for(set<BasicInteraction,BILess>::iterator iter_in=m_iset.begin();
	iter_in!=m_iset.end();
	iter_in++){
      Vec3 P1=m_bpart[iter_in->first()].getPos();
      Vec3 P2=m_bpart[iter_in->second()].getPos();
      if(iter->intersect(P1,P2)){ // on different sides
	set<BasicInteraction,BILess>::iterator h=iter_in;
	iter_in++;
	m_iset.erase(h);
	iter_in--;
      }
    }
  }
}

/*!
  Get closest line/line segment to a particle. Overloaded from  ARandomAssembly2D to include line segments with 
  overlap of 0.0

  \param Po the particle
*/
Line *FaultedBlock2D::getClosestPlane(const SimpleParticle& Po)
{
  //cout << "getClosestPlane : " << Po.getPos() << endl;
  Line* PL=ARandomAssembly2D::getClosestPlane(Po); // get closest boundary line

  Vec3 PoPos=Po.getPos();
  double dist=PL->sep(PoPos);
  //std::cout<< "posn: " << PoPos << std::endl;
  //cout << "plane: " << PL.GetO() << PL.GetN() << dist << endl;
  for(vector<pair<double,LineSegment> >::iterator iter=m_fault.begin();iter!=m_fault.end();iter++){
    double ndist=iter->second.sep(PoPos);
    //cout << "plane: " << (iter->second).GetO() << " - " << (iter->second).GetN() << " dist: " << ndist << endl;
    double dirdist=iter->second.GetN()*(PoPos-iter->second.GetO());
    if((ndist<dist) && (dirdist>0)){
      PL=&(iter->second);
      dist=ndist;
    }
  }
  //cout << " --- closest plane: " << PL->GetO() << " - "  << PL->GetN() << " dist: " << dist << endl;

  return PL;
}

/*!
  Tag particles along the split line

  \param tag1 the tag for particles "above" the split (y>y_split)
  \param tag2 the tag for particles "below" the split (y<y_split)
  \param d maximum distance from the split line at which a particle gets tagged
*/
void FaultedBlock2D::tagSplit(int tag1,int tag2,double d)
{
  cout << "FaultBlock2D::tagSplit" << endl;
  
  int ns=m_f2.size();
  for(int is=0;is<ns;is++){
    double rs=m_fault[is*2].first;
    for(vector<SimpleParticle>::iterator P_iter=m_bpart.begin();
	P_iter!=m_bpart.end();
	P_iter++){
      double di=m_f2[is].sep(P_iter->getPos());
      if(di<P_iter->getRad()+d+rs){
	double dir_sep=(P_iter->getPos()-m_f2[is].GetO())*m_f2[is].GetN();
	if(dir_sep>0.0){
	  P_iter->setTag(tag1);
	} else {
	  P_iter->setTag(tag2);
	}
      }
    }
   }
}
