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

#include "RandomBlock.h"

//-- IO includes --
#include <fstream>
#include <iomanip>
using std::ofstream;

//-- System includes --
#include <cstdlib>
using std::srand;

/*!
  Constructor of CRandomBlock

  \param xmin minimum in x-direction
  \param xmax maximum in x-direction
  \param ymin minimum in y-direction
  \param ymax maximum in y-direction
  \param rmin minimum particle radius
  \param rmax maximum particle radius
  \param mcd maximum relative distance for bond generation
*/
CRandomBlock2D::CRandomBlock2D(double xmin,double xmax,double ymin,double ymax,double rmin,double rmax,double mcd, bool circ_x)
{
  m_xmin=xmin;
  m_xmax=xmax;
  m_ymin=ymin;
  m_ymax=ymax;
  m_rmin=rmin;
  m_rmax=rmax;
  m_circ_x=circ_x;
  m_maxConnDist=mcd;
   // setup borders
  if(!circ_x) {
    Borders.push_back(Line(Vec3(1.0,0.0,0.0),Vec3(m_xmin,0.0,0.0)));  // left
    Borders.push_back(Line(Vec3(-1.0,0.0,0.0),Vec3(m_xmax,0.0,0.0))); // right
  }
  Borders.push_back(Line(Vec3(0.0,1.0,0.0),Vec3(0.0,m_ymin,0.0))); // bottom
  Borders.push_back(Line(Vec3(0.0,-1.0,0.0),Vec3(0.0,m_ymax,0.0))); // top
  // neighbor table
  m_snt=new CSimple2DNTable(Vec3(m_xmin,m_ymin,0.0),Vec3(m_xmax-m_xmin,m_ymax-m_ymin,0.0),2.1*rmax,circ_x);
}

CRandomBlock2D::~CRandomBlock2D()
{
  if(m_snt!=NULL) delete m_snt;
}

/*!
  Generate a random point within the space of the block
*/
Vec3 CRandomBlock2D::getAPoint()
{
  double px,py;
  
  px=m_random(m_xmin+m_rmin,m_xmax-m_rmin);
  py=m_random(m_ymin+m_rmin,m_ymax-m_rmin);

  return Vec3(px,py,0.0);
}

/*!
  Fill the space in the block

  \param tries number of times the insertion of a particle is tried 
  \param seed seed for the random number generator
*/
void CRandomBlock2D::generate(int tries,unsigned int seed)
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
  m_snt->getInteractions(m_iset,m_maxConnDist);
}

/*!
  Insert a particle into the internal structures

  \param P the particle
*/
void CRandomBlock2D::insertParticle(const SimpleParticle P)
{
  m_bpart.push_back(P);
  m_snt->insertParticle(P);
}

/*!
  Write the particles contained in the random block into a 
  LSM geometry file v 1.1

  \param filename the name of the file
*/
void CRandomBlock2D::writeToGeoFile(const string& filename)
{
  ofstream outfile;

  // open file
  outfile.open(filename.c_str());
  
  outfile << "LSMGeometry 1.2" << endl;

  // bounding box
  outfile << "BoundingBox " << m_xmin << " " << m_ymin << " -0.1 " << m_xmax << " " << m_ymax << " 0.1" << endl;
  outfile << "PeriodicBoundaries " << m_circ_x << " 0 0" << endl;
  outfile << "Dimension 2D" << endl;

  // particles
  outfile << "BeginParticles" << endl;
  outfile << "Simple" << endl;
  outfile << m_bpart.size() << endl;
	const int precision = 15;
  int count=0;
  for(vector<SimpleParticle>::const_iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
    outfile << std::setprecision(precision) << iter->getPos().X() << " " << iter->getPos().Y() << " " 
	    << iter->getPos().Z() << " " << iter->getRad() << " " 
	    << iter->getID() << " " << iter->getTag() << "\n";
    count++;
  }
  
  outfile << "EndParticles" << endl;
  // connections
  outfile << "BeginConnect" << endl;
  outfile << m_iset.size() << endl;
  for(set<BasicInteraction,BILess>::const_iterator iter=m_iset.begin();
      iter!=m_iset.end();
      iter++){
    outfile << *iter << "\n";
  }
  outfile << "EndConnect" << endl;
  // close file
  outfile.close();
}

/*!
  calculate the porosity of the material
*/
double CRandomBlock2D::calcPorosity()
{
  double a_total,a_part;
  
  a_total=(m_xmax-m_xmin)*(m_ymax-m_ymin);
  a_part=0.0;
  for(vector<SimpleParticle>::const_iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
    double radius=iter->getRad();
    a_part+=M_PI*radius*radius;
  }

  return 1.0-a_part/a_total;
}

/*!
  return a histogram of the particle size distribution

  \param nbins number of bins
*/
vector<pair<double,double> > CRandomBlock2D::getSizeDistribution(int nbins)
{
  vector<pair<double,double> >res(nbins);
  double inc=1.0/m_bpart.size();

  // setup bin centres
  for(int i=0;i<nbins;i++){
    double bc=m_rmin+(m_rmax-m_rmin)*((double(i)+0.5)/double(nbins));
    res[i].first=bc;
  }
  // calc distribution
  for(vector<SimpleParticle>::const_iterator iter=m_bpart.begin();
      iter!=m_bpart.end();
      iter++){
    double radius=iter->getRad();
    int bin=int(((radius-m_rmin)/(m_rmax+m_small_value-m_rmin))*double(nbins));
    res[bin].second+=inc;
  }

  return res;
}

/*!
  Tag particle closest to a given position

  \param pos the position
  \param tag the tag
*/
void CRandomBlock2D::tagParticleClosestTo(const Vec3& pos,int tag)
{
  int id=m_snt->getClosestParticleID(pos);
  vector<SimpleParticle>::iterator iter=m_bpart.begin();
  while((iter!=m_bpart.end()) && (iter->getID()!=id)) iter++;
  if(iter!=m_bpart.end()) iter->setTag(tag);
}

/*!
  Tag particles along x-edges

  \param tag1 tag for particles along x_min
  \param tag2 tag for particles along x_max
  \param d maximum distance from the edge at which a particle gets tagged
*/
void CRandomBlock2D::tagEdgeY(int tag1,int tag2,double d)
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
