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

#ifndef __ESYS_LSM_SPHAGGGOUGEBLOCK3D_H
#define __ESYS_LSM_SPHAGGGOUGEBLOCK3D_H

// --- project includes --
#include "Foundation/vec3.h"
#include "Geometry/GougeBlock3D.h"
#include "Geometry/SimpleParticle.h"
#include "Geometry/SphereBlockGenerator.h"

// --- STL includes ---
#include <vector>
using std::vector;

namespace esys {
  namespace lsm {

    /*!
      \class SphAggGougeBlock 
      \brief Block of gouge consisting of spherical aggregate grains  
    */
    class SphAggGougeBlock : public GougeBlock3D
    {
    public:
      typedef boost::shared_ptr<SphereBlockGenerator> SBG_ptr;

    protected:
      double m_min_rad_grain;
      double m_max_rad_grain;
      vector<SimpleParticle> m_macro_grains;
      NTablePtr m_nTablePtr2;
      ParticlePoolPtr m_particlePoolPtr2;
      GeneratorPtr m_grainGen;
      vector<SBG_ptr> m_grainParticleGen;
      int m_min_grain_tag;

      void generateMacroGrains();
      void fillMacroGrains();
      void setupNT2();
      void createInteractionSet();
      virtual void createGougeBlockGenerators();
      
    public:
      SphAggGougeBlock(const GougeBlockPrms&,double,double,int);
      virtual void generate();

      template <typename TmplVisitor> void visitParticles(TmplVisitor&);
      template <typename TmplVisitor> void visitParticles(TmplVisitor&) const;
    };

    /*!
      \class SphAggInteractionValidator
      \brief Used to check the validity of an interaction in a SphAggGougeBlock

      \author Steffen Abe
       $Revision$
       $Date$
    */ 
    class SphAggInteractionValidator
      {    
      private:
	const SphAggGougeBlock *m_pGougeBlock;
	double m_tolerance;
	int m_grain_tag;

      public:
	SphAggInteractionValidator(const SphAggGougeBlock&, double,int);
	bool isValid(const SimpleParticle&, const SimpleParticle&) const;
      };
  }
}

#include "SphAggGougeBlock.hpp"

#endif // __ESYS_LSM_SPHAGGGOUGEBLOCK3D_H 
