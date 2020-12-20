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

#include "Foundation/console.h"
#include "Geometry/RandomBlock3D.h"

//-- IO includes --
#include <fstream>
using std::ofstream;

//-- System includes --
#include <cstdlib>
using std::srand;

//-- project includes --
#include "SimpleNTable3D.h"

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
  \param mcd maximum relative distance for bond generation
*/
CRandomBlock3D::CRandomBlock3D(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double rmin,double rmax,double mcd,bool circ_x,bool is_bonded)
{
  m_xmin=xmin;
  m_xmax=xmax;
  m_ymin=ymin;
  m_ymax=ymax;
  m_zmin=zmin;
  m_zmax=zmax;
  m_rmin=rmin;
  m_rmax=rmax;
  m_circ_x=circ_x;
  m_is_bonded=is_bonded;
  m_maxConnDist=mcd;
   // setup borders
  if(!circ_x) {
    Borders.push_back(Plane3D(Vec3(1.0,0.0,0.0),Vec3(m_xmin,0.0,0.0)));  // left
    Borders.push_back(Plane3D(Vec3(-1.0,0.0,0.0),Vec3(m_xmax,0.0,0.0))); // right
  }
  Borders.push_back(Plane3D(Vec3(0.0,1.0,0.0),Vec3(0.0,m_ymin,0.0))); // bottom
  Borders.push_back(Plane3D(Vec3(0.0,-1.0,0.0),Vec3(0.0,m_ymax,0.0))); // top
  Borders.push_back(Plane3D(Vec3(0.0,0.0,1.0),Vec3(0.0,0.0,m_zmin))); // front
  Borders.push_back(Plane3D(Vec3(0.0,0.0,-1.0),Vec3(0.0,0.0,m_zmax))); // back
  // neighbor table
  m_snt=new CSimple3DNTable(Vec3(m_xmin,m_ymin,m_zmin),Vec3(m_xmax-m_xmin,m_ymax-m_ymin,m_zmax-m_zmin),2.1*rmax,circ_x);
}

CRandomBlock3D::~CRandomBlock3D()
{
  if(m_snt!=NULL) delete m_snt;
}


/*!
  Generate a random point within the space of the block
*/
Vec3 CRandomBlock3D::getAPoint()
{
  double px,py,pz;
  
  px=m_random(m_xmin+m_rmin,m_xmax-m_rmin);
  py=m_random(m_ymin+m_rmin,m_ymax-m_rmin);
  pz=m_random(m_zmin+m_rmin,m_zmax-m_rmin);

  return Vec3(px,py,pz);
}

/*!
  Fill the space in the block

  \param tries number of times the insertion of a particle is tried 
  \param seed seed for the random number generator
*/
void CRandomBlock3D::generate(int tries,unsigned int seed)
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
  if(m_is_bonded){
    m_snt->getInteractions(m_iset,m_maxConnDist);
  }
}

/*!
  Insert a particle into the internal structures

  \param P the particle
*/
void CRandomBlock3D::insertParticle(const SimpleParticle P)
{
  m_bpart.push_back(P);
  m_snt->insertParticle(P);
}

/*!
  Tag particle closest to a given position

  \param pos the position
  \param tag the tag
*/
//void CRandomBlock3D::tagParticleClosestTo(const Vec3& pos,int tag)
//{}

/*!
  Tag particles along xz-edges

  \param tag1 tag for particles along y_min
  \param tag2 tag for particles along y_max
  \param d maximum distance from the edge at which a particle gets tagged
*/
void CRandomBlock3D::tagEdgeY(int tag1,int tag2,double d)
{
  for(vector<SimpleParticle>::iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
    double py=iter->getPos().Y();
    double r=iter->getRad();
    if(py-m_ymin<r+d) iter->setTag(tag1);
    if(m_ymax-py<r+d) iter->setTag(tag2);
  }
}
 
/*!
  Tag particles along xy-edges

  \param tag1 tag for particles along z_min
  \param tag2 tag for particles along z_max
  \param d maximum distance from the edge at which a particle gets tagged
*/
void CRandomBlock3D::tagEdgeZ(int tag1,int tag2,double d)
{
  for(vector<SimpleParticle>::iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
    double py=iter->getPos().Z();
    double r=iter->getRad();
    if(py-m_zmin<r+d) iter->setTag(tag1);
    if(m_zmax-py<r+d) iter->setTag(tag2);
  }
}
 
/*!
  Write the particles contained in the random block into a 
  LSM geometry file v 1.2

  \param filename the name of the file
*/
void CRandomBlock3D::writeToGeoFile(const string& filename)
{
  ofstream outfile;

  // open file
  outfile.open(filename.c_str());
  
  outfile.precision(10);

  outfile << "LSMGeometry 1.2" << endl;

  // bounding box
  outfile << "BoundingBox " << m_xmin << " " << m_ymin << " " << m_zmin << " " 
	  << m_xmax << " " << m_ymax << " " << m_zmax << endl;
  outfile << "PeriodicBoundaries " << m_circ_x << " 0 0" << endl;
  outfile << "Dimension 3D" << endl;
 
  // particles
  outfile << "BeginParticles" << endl;
  outfile << "Simple" << endl;
  outfile << m_bpart.size() << endl;
  int count=0;
  for(vector<SimpleParticle>::const_iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
    outfile << iter->getPos().X() << " " << iter->getPos().Y() << " " 
	    << iter->getPos().Z() << " " << iter->getRad() << " " 
	    << iter->getID() << " " << iter->getTag() << endl;
    count++;
  }
  
  outfile << "EndParticles" << endl;
  // connections
  outfile << "BeginConnect" << endl;
  outfile << m_iset.size() << endl;
  for(set<BasicInteraction,BILess>::const_iterator iter=m_iset.begin();
      iter!=m_iset.end();
      iter++){
    outfile << *iter << endl;
  }
  outfile << "EndConnect" << endl;
  // close file
  outfile.close();
}

/*!
  calculate the porosity of the material
*/
double CRandomBlock3D::calcPorosity()
{
  double v_total,v_part;

  v_total=(m_xmax-m_xmin)*(m_ymax-m_ymin)*(m_zmax-m_zmin);
  v_part=0.0;
  for(vector<SimpleParticle>::const_iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
    double radius=iter->getRad();
    v_part+=(4.0/3.0)*M_PI*radius*radius*radius;
  }
  
  console.Info() << "total volume : " << v_total << "\n";
  console.Info() << "filled volume: " << v_part << "\n";
  return 1.0-v_part/v_total;
}

/*!
  return a histogram of the particle size distribution

  \param nbins number of bins
*/
vector<pair<double,double> > CRandomBlock3D::getSizeDistribution(int nbins)
{
  vector<pair<double,double> >res(nbins);

  return res;
}
