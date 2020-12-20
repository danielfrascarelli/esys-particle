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

#include "Geometry/SphereBlockGenerator.h"

// -- project includes ---
#include "Foundation/BoundingBox.h"
#include "Geometry/GridIterator.h"
#include "Geometry/Sphere3d.h"

// -- system include ---
#include <cstdlib>
using std::rand;

namespace esys
{
  namespace lsm
  {  
    /*!
      constructor

      \param ntable the neigbour table to be used 
      \param pool the particle pool
      \param tol the fitting tolerance
      \param pos center position
      \param rad radius
      \param rmin minimum particle radius
      \param rmax maximum particle radius
      \param ntries max. nr. of tries
    */
    SphereBlockGenerator::SphereBlockGenerator(NTable& ntable,ParticlePool& pool,double tol,const Vec3& pos,double rad,double rmin,double rmax,int ntries,int tag)
      :ParticleGenerator(ntable, pool)
    {
      m_tol=tol;
      m_center=pos;
      m_radius=rad;
      m_min_rad=rmin;
      m_max_rad=rmax;
      m_max_tries=ntries;
      m_tag=tag;
    }

    /*!
      destructor
    */
    SphereBlockGenerator::~SphereBlockGenerator()
    {}

    /*!
      get next available particle ID
    */
    int SphereBlockGenerator::getNextId()
    {
      return static_cast<int>(getNTable().getNumParticles());
    }
    

    /*!
      generate a particle at a given position

      \param point
    */
    SimpleParticle SphereBlockGenerator::generateParticle(const Vec3 &point)
    {
      return SimpleParticle(point, getRadius(), getNextId(),m_tag);
    }

    /*!
      get a random point inside the sphere
    */
    Vec3 SphereBlockGenerator::getAPoint()
    {
      Vec3 res;

      do{
	double x=m_radius*(1.0-2.0*((static_cast<double>(rand()))/(static_cast<double>(RAND_MAX))));
	double y=m_radius*(1.0-2.0*((static_cast<double>(rand()))/(static_cast<double>(RAND_MAX))));
	double z=m_radius*(1.0-2.0*((static_cast<double>(rand()))/(static_cast<double>(RAND_MAX))));
	res=Vec3(x,y,z);
      } while (res.norm()>=m_radius);

      return m_center+res;
    }

    /*!
      ??
    */
    double SphereBlockGenerator::getRadius() const
    {
      return m_min_rad+((m_max_rad-m_min_rad)*(static_cast<double>(rand()))/(static_cast<double>(RAND_MAX)));
    }

    /*!
      return max. radius to be used as spacing for grid iterator
    */
    double SphereBlockGenerator::getGridRadius() const
    {
      return m_max_rad;
    }

    /*!
      check particle fit
    */
    bool SphereBlockGenerator::particleFits(const SimpleParticle &particle) const
    {
      bool radius_fit=(particle.getRad() >= m_min_rad) && (particle.getRad() <= m_max_rad);
      bool inside=((m_center-particle.getPos()).norm()<=(m_radius-particle.getRad()));

      return radius_fit && inside;
    }

    /*! 
      calculate & return bounding box for bounding sphere
    */
    const BoundingBox SphereBlockGenerator::getBBox() const
    {
      Vec3 minvec=m_center-Vec3(m_radius,m_radius,m_radius);
      Vec3 maxvec=m_center+Vec3(m_radius,m_radius,m_radius);
      return BoundingBox(minvec,maxvec);
    }
    
    /*!
      insert particle 

      \param particle
    */
    void SphereBlockGenerator::insertParticle(const SimpleParticle &particle)
    {
      SimpleParticle *pParticle = getParticlePool().construct(particle);
      pParticle->setTag(m_tag);
      m_particleVector.push_back(pParticle);
      m_idSet.insert(pParticle->getID());
      getNTable().insert(pParticle);
    }
    /*!
      generate particles
    */
    void SphereBlockGenerator::generate()
    {
      generateSeedParticles();
      generateFillParticles();
    }

    /*!
      generate seed particles
    */
    void SphereBlockGenerator::generateSeedParticles()
    {
      // setup grid iterator
      GridIterator pointIt = GridIterator(getBBox(), getGridRadius());
      // run through grid
      while (pointIt.hasNext()) {
        SimpleParticle particle = generateParticle(pointIt.next());
        if (particleFits(particle)) {
          insertParticle(particle);
        }
      }
    }
    
    /*!
      get closes Neigbours
      
      \param P the particle
      \param np the number of neighbours
    */
    vector<SimpleParticle*> SphereBlockGenerator::getClosestNeighbors(const SimpleParticle& P,int np)
    {
      vector<SimpleParticle*> nv=getNTable().getUniqueNeighbourVector(P.getPos(), m_max_rad + m_tol);
      std::sort(nv.begin(), nv.end(), ParticleComparer(P));
      if (nv.size() > static_cast<size_t>(np)) {
        nv.erase(nv.begin() + np, nv.end());
      }
      return nv;
    }

    /*!
      Find a fit for a sphere using the list of neigbors and the outer sphere
      
      \param Po the particle to fit
      \param NL the list of neighbors
    */
    bool SphereBlockGenerator::findAFitWithSphere(SimpleParticle& Po, const vector<SimpleParticle*>& NL)
    {
      //      cout << "findAFit - 3 particles + outer\n";
      bool find_a_fit ;
      Vec3 M;
      double r;
      int id=Po.getID();

      if(NL.size()<3){
	find_a_fit=false;
	//	cout << "less than 3 neighbors" << endl; // can't happen
      } else {
	Vec3 Pos1=m_center;
	Vec3 Pos2=NL[0]->getPos();
	Vec3 Pos3=NL[1]->getPos();
	Vec3 Pos4=NL[2]->getPos();
	double r1=-1.0*m_radius;
	double r2=NL[0]->getRad();
	double r3=NL[1]->getRad();
	double r4=NL[2]->getRad();

	find_a_fit=Sphere3D::FillIn(Pos1,Pos2,Pos3,Pos4,r1,r2,r3,r4,M,r);
	Po=SimpleParticle(M,r,id,m_tag);
	//	cout << "found " << M << " , " << r << endl; 
      }
   
      return find_a_fit ;
    }

    /*!
      Find a fit for a sphere using the list of neigbors
      
      \param Po the particle to fit
      \param NL the list of neighbors
    */
    bool SphereBlockGenerator::findAFit(SimpleParticle& Po, const vector<SimpleParticle*>& NL)
    {
//       cout << "findAFit - 4 particles\n";
      bool find_a_fit ;
      Vec3 M;
      double r;
      int id=Po.getID();

      if(NL.size()<4){
	find_a_fit=false;
// 	cout << "less than 4 neighbors" << endl; // can't happen
      } else {
	Vec3 Pos1=NL[0]->getPos();
	Vec3 Pos2=NL[1]->getPos();
	Vec3 Pos3=NL[2]->getPos();
	Vec3 Pos4=NL[3]->getPos();
	double r1=NL[0]->getRad();
	double r2=NL[1]->getRad();
	double r3=NL[2]->getRad();
	double r4=NL[3]->getRad();

	find_a_fit=Sphere3D::FillIn(Pos1,Pos2,Pos3,Pos4,r1,r2,r3,r4,M,r);
	Po=SimpleParticle(M,r,id,m_tag);
// 	cout << "found " << M << " , " << r << endl; 
      }
   
      return find_a_fit ;
    }

    /*!
      check if Po is within the Space and is not crossing any boundary or 
      overlapping with other particles.
      
      \param Po the particle
    */
    bool SphereBlockGenerator::checkAFit(const SimpleParticle& Po)
    {
      bool fail=false;
	
      // check vs. radius
      fail=!particleFits(Po);
      // check vs. all neighbors
      if(!fail){
	vector<SimpleParticle*> NL=getNTable().getUniqueNeighbourVector(Po.getPos(), m_max_rad + m_tol);; 
	vector<SimpleParticle*>::const_iterator iter=NL.begin();
	while(!fail && iter!=NL.end()){
	  double dist=(Po.getPos()-(*iter)->getPos()).norm()+m_tol;
	  if(dist < (Po.getRad()+(*iter)->getRad())){
	    fail=true;
	    //cout << "Fail : particle collision" << endl;
	  }
	  iter++;
	}
      }
	
      return !fail ;
    }
      
    /*!
      fill in 
    */
    void SphereBlockGenerator::generateFillParticles()
    {
      int countfail=0;
      int countfound=0;
      //bool fail,findfit;
      bool findfit=false;
      //bool foundwithsphere=false;

      cout << "SphereBlockGenerator::generateFillParticles" << endl;
      while(countfail<m_max_tries){
	findfit=false;
	Vec3 P=getAPoint();
	//	cout << "Got Point: " << P << endl;
	SimpleParticle Po=generateParticle(P);
	vector<SimpleParticle*> T4=getClosestNeighbors(Po,4);
	if(T4.size()>3){ // at least 4 neighbors 
	  SimpleParticle* Pi=T4[0];
	  double ndist=(Po.getPos()-Pi->getPos()).norm();
	  if( ndist==0.0){
	    findfit=false;
	  } else {
	    if( ndist < Pi->getRad()){ // if Po inside Pi -> move Po to the surface of Pi
	      Vec3 npos=Pi->getPos()+((Po.getPos()-Pi->getPos())*(Pi->getRad()/ndist));
	      Po.moveTo(npos);
	    }
	    double dist_p=m_radius-(Po.getRad()+(m_center-Po.getPos()).norm()); // distance between Po and outer sphere
	    double dist_3=(Po.getPos()-T4[3]->getPos()).norm()-T4[3]->getRad();
	    if (dist_p>dist_3){  // 4th particle closer than plane -> fit 4 particles
	      findfit=findAFit(Po,T4);
	    } else { // outer sphere closer than 4th particle -> fit 3 particles + outer sphere
	      findfit=findAFitWithSphere(Po,T4);
	    }
	  }
	} else {  // 3 neighbors  -> try 3 particles + outer sphere
	  SimpleParticle* Pi=T4[0];
	  double ndist=(Po.getPos()-Pi->getPos()).norm();
	  if( ndist==0.0){
	    findfit=false;
	  } else {
	    if( ndist < Pi->getRad()){ // if Po inside Pi -> move Po to the surface of Pi
	      Vec3 npos=Pi->getPos()+((Po.getPos()-Pi->getPos())*(Pi->getRad()/ndist));
	      Po.moveTo(npos);
	    }
	    findfit=findAFitWithSphere(Po,T4);
	  }
	}
	if(findfit){ // found something, check
	  findfit=checkAFit(Po);
	}
	if(findfit){  // found & checked -> insert
	  insertParticle(Po);
	  if(countfail*10>m_max_tries){
	    cout << "found particle " << Po.getID() << " after " << countfail << " tries" << endl;
	  }
	  countfail=0;
	  countfound++;
	      } else {
	  countfail++;
	}
      }
      cout << "end SphereBlockGenerator::generateFillParticles" << endl;
    }

  }
}
