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

#include "Geometry/GranularGougeBlock3D.h"
#include "Foundation/BoundingBox.h"

// -- I/O includes --
#include <iostream>

namespace esys 
{
  namespace lsm {
    
    //======== GranularGougeBlock3D ==============

    /*!
      Constructor for GranularGougeBlock3D. Do nothing and call the base class
      constructor (GougeBlock3D)

      \param prms the GougeBlock3D parameters
    */
    GranularGougeBlock3D::GranularGougeBlock3D(const GougeBlockPrms &prms) 
      :GougeBlock3D(prms)
    {
    }

    /*!
      Destructor. No dynamically allocated data in class -> do nothing 
    */
    GranularGougeBlock3D::~GranularGougeBlock3D()
    {
    }

    /*!
      Create the seed points for the grain generation algorithm. Algorithm currently assumes single
      gouge layer.

      \param sdx seed density, i.e. average distance between seeds in x-direction
      \param sdy seed density in y-direction
      \param sdz seed density in z-direction
      \param rdx random variation of seed points in x-direction
      \param rdy random variation of seed points in y-direction
      \param rdz random variation of seed points in z-direction
    */
    void GranularGougeBlock3D::generateSeeds(double sdx,double sdy,double sdz,double rdx,double rdy,double rdz)
    {
      // get bouding box
      PackingInfoVector packingInfoVector = m_prms.getGougePackingInfoVector();
      BoundingBox bbx=(packingInfoVector.begin())->getBBox();
      Vec3 bbx_min_pt=bbx.getMinPt();
      Vec3 bbx_dims=bbx.getMaxPt()-bbx_min_pt;

      // calculate no. of seeds per direction
      int n_seed_x=int(floor(bbx_dims.X()/sdx));
      int n_seed_y=int(floor(bbx_dims.Y()/sdy));
      int n_seed_z=int(floor(bbx_dims.Z()/sdz));

      // generate them
      for(int i=0;i<n_seed_x;i++){
	for(int j=0;j<n_seed_y;j++){
	  for(int k=0;k<n_seed_z;k++){
	    Vec3 seedpos=Vec3((double(i)+0.5)*sdx,(double(j)+0.5)*sdy,(double(k)+0.5)*sdz);
	    double rx=rdx*(-0.5+((double)(rand())/(double)(RAND_MAX)));
	    double ry=rdy*(-0.5+((double)(rand())/(double)(RAND_MAX)));
	    double rz=rdz*(-0.5+((double)(rand())/(double)(RAND_MAX)));
	    m_grain_seeds.push_back(bbx_min_pt+seedpos+Vec3(rx,ry,rz));
	  }
	}
      }
    }

    /*!
      Generate composite grains from the existing gouge particles by randomly distributing
      seed points within the gouge region, then tagging all particles closest to the same 
      seed point with the same tag.

      \warning currently assumes single gouge layer/generator

      \param sdx seed density, i.e. average distance between seeds in x-direction
      \param sdy seed density in y-direction
      \param sdz seed density in z-direction
      \param rdx random variation of seed points in x-direction
      \param rdy random variation of seed points in y-direction
      \param rdz random variation of seed points in z-direction
      \param min_tag minimum tag to be used in order not to collide with allready used tags
      \param rm_threshold grains with less then rm_threshold particles get removed. Defaults to 0
    */
    //void GranularGougeBlock3D::generateGrains(double sdx,double sdy,double sdz,double rdx,double rdy,double rdz,int min_tag,int rm_threshold)
    void GranularGougeBlock3D::generateGrains(double sdx,double sdy,double sdz,double rdx,double rdy,double rdz,int min_tag)
    {
      generateSeeds(sdx,sdy,sdz,rdx,rdy,rdz);
      for(vector<Vec3>::iterator iter=m_grain_seeds.begin();
	  iter!=m_grain_seeds.end();
	  iter++){
	std::cout << *iter << std::endl;
      }
      // make grains
      GeneratorPtrVector::iterator it = m_gougeGenPtrVector.begin();
      BlockGenerator::ParticleIterator particleIt = (*it)->getParticleIterator();
      int n_seeds=m_grain_seeds.size();
      while (particleIt.hasNext()) {
	SimpleParticle* current_particle=particleIt.next();
	double nearest_seed_dist=(current_particle->getPos()-m_grain_seeds[0]).norm();
	int nearest_seed_id=0;
	for(int i=1;i<n_seeds;i++){ // assume there are at least 2 seeds
	  double dist=(current_particle->getPos()-m_grain_seeds[i]).norm();
	  if(dist<nearest_seed_dist){
	    nearest_seed_dist=dist;
	    nearest_seed_id=i;
	  }
	}
	current_particle->setTag(nearest_seed_id+min_tag);
      }
    }
	
    /*!
      Create interaction set. Changed from base class by using a different validator which
      allows links between particles with the same tag, i.e. belonging to the same composite grain.
      Refactor ?
    */ 
    void GranularGougeBlock3D::createInteractionSet()
    {
      NTable::ParticleIterator particleIt = m_nTablePtr->getParticleIterator();
      GranularInteractionValidator validator = GranularInteractionValidator(*this, m_prms.getConnectionTolerance());
      while (particleIt.hasNext()) {
        const NTable::Particle *pParticle = particleIt.next();
        const NTable::ParticleVector neighbours =
          m_nTablePtr->getNeighbourVector(
            pParticle->getPos(),
            pParticle->getRad() + m_prms.getConnectionTolerance()
          );
        for ( NTable::ParticleVector::const_iterator it = neighbours.begin();
	      it != neighbours.end();
	      it++ ){
          if (validator.isValid(*pParticle, *(*it))) {
	    //	    cout << "inserting " << pParticle->getID() << "-" << (*it)->getID() << endl;; 
            m_interactionSet.insert(BasicInteraction(pParticle->getID(), (*it)->getID()));
          }
        }
      }
    }

    /*!
     */
    void GranularGougeBlock3D::generate()
    {
      for (
        GeneratorPtrVector::iterator it = m_genPtrVector.begin();
        it != m_genPtrVector.end();
        it++
      )
      {
        (*it)->generate();
      }
    }
 
    //========= GranularInteractionValidator ============
    
    /*!
     */
    GranularInteractionValidator::GranularInteractionValidator(const GranularGougeBlock3D& gougeBlock, double tolerance)
      : m_pGougeBlock(&gougeBlock),
	m_tolerance(tolerance)
    {}
    
    /*!
     */
    bool GranularInteractionValidator::isValid(const SimpleParticle &p1, const SimpleParticle &p2) const
    {
       return (
	      (p1.getID() < p2.getID()) // order correct -> avoid doubles
	      &&
	      ((p1.getPos() - p2.getPos()).norm() < (m_tolerance + (p1.getRad() + p2.getRad()))) // distance correct
	      &&
	      ((((!m_pGougeBlock->isGougeParticle(p1)) && (!m_pGougeBlock->isGougeParticle(p2)))
	       &&
	       ((!m_pGougeBlock->areInDifferentFaultBlocks(p1, p2)))) // both outside an on the same side of the gouge
	      ||
	       (((m_pGougeBlock->isGougeParticle(p1)) && (m_pGougeBlock->isGougeParticle(p2)))
	       &&
	      (p1.getTag()==p2.getTag()))) // both inside and same tag, i.e. same grain
	      );
    }
  }
}
