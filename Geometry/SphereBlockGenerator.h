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

#ifndef __ESYS_LSM_SPHEREBLOCKGENERATOR
#define __ESYS_LSM_SPHEREBLOCKGENERATOR

// --- project includes ---
#include <Geometry/ParticleGenerator.h>
#include <Geometry/SimpleParticle.h>

// --- STL includes ---
#include <set>
using std::set;

namespace esys
{
  namespace lsm
  {    
    /*!
      \class SphereBlockGenerator
    */
    class SphereBlockGenerator : public ParticleGenerator
    {
    public:
      // types 
      typedef NTable::ParticleVector ParticleVector;
      typedef NTable::ParticleIterator ParticleIterator;
      typedef set<int> IdSet;

      // functions
      SphereBlockGenerator(NTable&,ParticlePool&,double,const Vec3&,double,double,double,int,int);
      virtual ~SphereBlockGenerator();

      virtual void generate();
      virtual void generateSeedParticles();
      virtual void generateFillParticles();
      virtual SimpleParticle generateParticle(const Vec3 &point);
      virtual void insertParticle(const SimpleParticle&);
      virtual double getRadius() const;
      int getNextId();
      size_t getNumParticles() const {return m_idSet.size();};
      const BoundingBox getBBox() const;
      virtual double getGridRadius() const;
      virtual bool particleFits(const SimpleParticle &particle) const;
      ParticleIterator getParticleIterator(){return ParticleIterator(m_particleVector);}
      vector<SimpleParticle*> getClosestNeighbors(const SimpleParticle&,int);
      bool findAFitWithSphere(SimpleParticle&, const vector<SimpleParticle*>&);
      bool findAFit(SimpleParticle&, const vector<SimpleParticle*>&);
      bool checkAFit(const SimpleParticle&);
      Vec3 getAPoint();

    private:
      ParticleVector m_particleVector;
      double         m_tol;
      IdSet          m_idSet;
      Vec3           m_center;
      double         m_radius;
      double         m_min_rad;
      double         m_max_rad;
      int            m_max_tries;
      int            m_tag;
    };
  }
}

#endif // __ESYS_LSM_SPHEREBLOCKGENERATOR
