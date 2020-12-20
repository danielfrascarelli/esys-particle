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

#include "Geometry/SphAggGougeBlock.h"
#include "Geometry/RandomBlockGenerator.h"

// --- IO includes ---
#include <iostream>
using std::cout;
using std::endl;

namespace esys
{
  namespace lsm
  {
    /*!
      constructor

      \param params the general gouge parameters
      \param minRadGrain minimum radius of the macro grains
      \param maxRadGrain maximum radius of the macro grains
      \param minGrainTag minimum tag for grains 
    */
    SphAggGougeBlock::SphAggGougeBlock(const GougeBlockPrms& params,double minRadGrain,double maxRadGrain,int minGrainTag):
      GougeBlock3D(params),
      m_particlePoolPtr2(new ParticlePool(2048))
    {
      m_min_rad_grain=minRadGrain;
      m_max_rad_grain=maxRadGrain;
      m_min_grain_tag=minGrainTag;
    }

    /*!
      helper function for generate - setup inital sphere packing for macro grains
    */
    void SphAggGougeBlock::generateMacroGrains()
    {
      cout << "begin SphAggGougeBlock::generateMacroGrains" << endl;
      
      // get fault packing info - should check if it is actually there (size()==2)
      PackingInfoVector packingInfoVector = m_prms.getFaultPackingInfoVector();

      // get the bounding boxes of the fault blocks 
      BoundingBox bbx1=packingInfoVector[0].getBBox();
      BoundingBox bbx2=packingInfoVector[1].getBBox();

      // get corners of bounding boxes
      Vec3 minPt1=bbx1.getMinPt();
      Vec3 maxPt1=bbx1.getMaxPt();
      Vec3 minPt2=bbx2.getMinPt();
      Vec3 maxPt2=bbx2.getMaxPt();

      // check which points are the "inner" points
      double diff1=fabs(minPt1.Y()-maxPt2.Y());
      double diff2=fabs(minPt2.Y()-maxPt1.Y());
      double y1,y2;
      if(diff1<diff2){
	y1=min( minPt1.Y(),maxPt2.Y());
	y2=max( minPt1.Y(),maxPt2.Y());
      } else {
	y1=min( minPt2.Y(),maxPt1.Y());
	y2=max( minPt2.Y(),maxPt1.Y());
      }
      
      // generate bounding box for gouge
      Vec3 minPt=Vec3(minPt1.X(),y1,minPt1.Z());
      Vec3 maxPt=Vec3(maxPt1.X(),y2,maxPt1.Z());
       BoundingBox bbx3=BoundingBox(minPt,maxPt);
       cout << "bounding box for gouge: " << bbx3 << endl;

       // generate vector of bounding planes
       PlaneVector pv;
       pv.push_back(Plane3D(Vec3(0, 0,  1), bbx3.getMinPt()));
       pv.push_back(Plane3D(Vec3(0, 1,  0), bbx3.getMinPt()));
       pv.push_back(Plane3D(Vec3(1, 0,  0), bbx3.getMinPt()));
        pv.push_back(Plane3D(Vec3(0, 0,  -1), bbx3.getMaxPt()));
       pv.push_back(Plane3D(Vec3(0, -1,  0), bbx3.getMaxPt()));
       pv.push_back(Plane3D(Vec3(-1, 0,  0), bbx3.getMaxPt()));

       // setup grain generator
       m_grainGen =  
	 GeneratorPtr(
		      new RandomBlockGenerator(*(m_nTablePtr2.get()),
					       *(m_particlePoolPtr2.get()),
					       bbx3,
					       BoolVector(3,false),
					       m_prms.getTolerance(),
					       m_min_rad_grain,
					       m_max_rad_grain,
					       pv,
					       m_prms.getMaxInsertionFailures()
					       )
		      );
       m_grainGen->generate();
       cout << m_grainGen->getNumParticles() << " generated" << endl;
       
       cout << "end SphAggGougeBlock::generateMacroGrains" << endl;
    }

    /*!
      helper function for generate - fill each generated macro sphere with smaller particles
    */
    void SphAggGougeBlock::fillMacroGrains()
    {
      cout << "begin SphAggGougeBlock::fillMacroGrains" << endl;
      // --- setup a generator per grain
      // get min/max particle radius ( single gouge layer assumed )
      PackingInfoVector::const_iterator it = (m_prms.getGougePackingInfoVector()).begin();
      double rmin=it->getMinRadius();
      double rmax=it->getMaxRadius();
      // for all grains ....
      NTable::ParticleIterator particleIt = m_nTablePtr2->getParticleIterator();
      int curr_tag=m_min_grain_tag;
      while (particleIt.hasNext()) {
	SimpleParticle* SP=particleIt.next();
        cout << *SP << endl;
	SBG_ptr genPtr = SBG_ptr ( new SphereBlockGenerator (*(m_nTablePtr.get()),
							     (*(m_particlePoolPtr.get())),
							     m_prms.getTolerance(),
							     SP->getPos(),
							     SP->getRad(),
							     rmin,
							     rmax,
							     m_prms.getMaxInsertionFailures(),
							     curr_tag));
	m_grainParticleGen.push_back(genPtr);
	curr_tag++;
      }
      // use generators
      curr_tag=m_min_grain_tag;
      for (vector<SBG_ptr>::iterator it = m_grainParticleGen.begin();
	   it != m_grainParticleGen.end();
	   it++)
      {
	// generate particles
        (*it)->generate();
	// tag them
	SphereBlockGenerator::ParticleIterator particleIt = (*it)->getParticleIterator();
	while (particleIt.hasNext()) {
	  particleIt.next()->setTag(curr_tag);
	}
	cout << (*it)->getNumParticles() << "with tag " << curr_tag << " generated " << endl;
	curr_tag++;
      }
      cout << "end SphAggGougeBlock::fillMacroGrains" << endl;
    }

    /*!
      helper function for generate - setup the 2nd NTable (for the macro grains)
    */ 
    void SphAggGougeBlock::setupNT2()
    {
      const BoundingBox bBox = m_prms.getBBox();
      Vec3 ntableAdjust = Vec3(0.0,0.0,0.0);
      if (m_prms.getBBox().getSizes().Z() >= 4*m_prms.getMaxRadius()) {
        ntableAdjust += Vec3(m_prms.getMaxRadius(), 0, 0);
      }
      const BoundingBox nTableBBox(bBox.getMinPt(), bBox.getMaxPt() - ntableAdjust);
      m_nTablePtr2 =
        NTablePtr(
          new NTable(
            nTableBBox,
            (4.0*m_min_rad_grain), // grid spacing
            vector<bool>(3,false),   // non periodic bcond
            2.1*m_max_rad_grain // width of border-region in which 
                                // particles are duplicated
                                // for circular boundry
	    )
	  );
    }
    
    /*!
      create block generators for gouge region -> do nothing
    */
    void SphAggGougeBlock::createGougeBlockGenerators()
    {}
    
    /*!
      Create interaction set. Changed from base class by using a different validator which
      allows links between particles with the same tag, i.e. belonging to the same composite grain.
      Refactor ?
    */ 
    void SphAggGougeBlock::createInteractionSet()
    {
      NTable::ParticleIterator particleIt = m_nTablePtr->getParticleIterator();
      SphAggInteractionValidator validator = SphAggInteractionValidator(*this, m_prms.getConnectionTolerance(),m_min_grain_tag);
      while (particleIt.hasNext()) {
        const NTable::Particle *pParticle = particleIt.next();
        const NTable::ParticleVector neighbours =
          m_nTablePtr->getUniqueNeighbourVector(
            pParticle->getPos(),
            pParticle->getRad() + m_prms.getConnectionTolerance()
          );
        for ( NTable::ParticleVector::const_iterator it = neighbours.begin();
	      it != neighbours.end();
	      it++ ){
          if (validator.isValid(*pParticle, *(*it))) {
            m_interactionSet.insert(BasicInteraction(pParticle->getID(), (*it)->getID()));
          }
        }
      }
    }
    /*!
      generate particle packing
    */
    void SphAggGougeBlock::generate()
    {
      cout << "begin SphAggGougeBlock::generate" << endl;
      // --- generate solid parts ---
      // setup block generators
      createRegularBlockGenerators();
      createFaultBlockGenerators();
      
      // use block generators  
      std::cout << "bbox = " << m_prms.getBBox() << std::endl;
      for (
        GeneratorPtrVector::iterator it = m_genPtrVector.begin();
        it != m_genPtrVector.end();
        it++
      )
      {
        (*it)->generate();
      }
      // --- generate gouge ---
      setupNT2();
      generateMacroGrains();
      fillMacroGrains();
      // --- generate interactions ---
      createInteractionSet();
	    
      cout << "end SphAggGougeBlock::generate" << endl;
    }

    //========= GranularInteractionValidator ============
    
    /*!
     */
    SphAggInteractionValidator::SphAggInteractionValidator(const SphAggGougeBlock& gougeBlock, double tolerance, int graintag)
      : m_pGougeBlock(&gougeBlock),
	m_tolerance(tolerance),
	m_grain_tag(graintag)
    {}
    
    /*!
     */
    bool SphAggInteractionValidator::isValid(const SimpleParticle &p1, const SimpleParticle &p2) const
    {
      bool is_in_order=(p1.getID() < p2.getID());
      bool is_distance=((p1.getPos() - p2.getPos()).norm() < (m_tolerance + (p1.getRad() + p2.getRad())));
      bool is_same_side=(p1.getTag()<m_grain_tag) && (p2.getTag()<m_grain_tag);// && (!m_pGougeBlock->areInDifferentFaultBlocks(p1, p2));
      bool is_same_grain=(p1.getTag()>=m_grain_tag) && (p2.getTag()>=m_grain_tag) && (p1.getTag()==p2.getTag());
      // bool both_grain=(p1.getTag()>=m_grain_tag) && (p2.getTag()>=m_grain_tag);

      return (is_in_order && is_distance) && (is_same_side || is_same_grain);
    }
  }
}
