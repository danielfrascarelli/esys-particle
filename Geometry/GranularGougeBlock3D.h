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

#ifndef __ESYS_LSM_GRANULARGOUGEBLOCK3D_H
#define __ESYS_LSM_GRANULARGOUGEBLOCK3D_H

// --- project includes --
#include "Geometry/GougeBlock3D.h"

namespace esys {
  namespace lsm {

    /*!
      \class GranularGougeBlock3D
      \brief Class to generate a 3d block of material consisting of
      a granular gouge between two solid blocks

      \author Steffen Abe
       $Revision$
       $Date$
    */
    class GranularGougeBlock3D : public GougeBlock3D
      {
      private:
	vector<Vec3> m_grain_seeds; //!< seed points for grain generation algorithm

	void generateSeeds(double,double,double,double,double,double); 

      public:
	GranularGougeBlock3D(const GougeBlockPrms &prms);
	virtual ~GranularGougeBlock3D();

	virtual void createInteractionSet();
	virtual void generate();
	//virtual void generateGrains(double,double,double,double,double,double,int,int rm_threshold=0);
	virtual void generateGrains(double,double,double,double,double,double,int);
      };

    /*!
      \class GranularInteractionValidator
      \brief Used to check the validity of an interaction in a GranularGougeBlock

      \author Steffen Abe
       $Revision$
       $Date$
    */ 
    class GranularInteractionValidator
      {    
      private:
	const GranularGougeBlock3D *m_pGougeBlock;
	double m_tolerance;

      public:
	GranularInteractionValidator(const GranularGougeBlock3D&, double);
	bool isValid(const SimpleParticle&, const SimpleParticle&) const;
      };
  }
}

#endif // __ESYS_LSM_GRANULARGOUGEBLOCK3D_H
