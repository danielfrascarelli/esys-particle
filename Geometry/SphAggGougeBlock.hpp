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

namespace esys {
  namespace lsm {

    /*!
     */
    template <typename TmplVisitor> 
    void SphAggGougeBlock::visitParticles(TmplVisitor& visitor)
    {
      GougeBlock3D::visitParticles(visitor);
      for(vector<SBG_ptr>::iterator it = m_grainParticleGen.begin();
	  it != m_grainParticleGen.end();
	  it++)
	{
	  SphereBlockGenerator::ParticleIterator particleIt = (*it)->getParticleIterator();
          while (particleIt.hasNext()) {
            particleIt.next()->visit(visitor);
          }
	}
    }

    /*!
     */
    template <typename TmplVisitor> 
    void SphAggGougeBlock::visitParticles(TmplVisitor& visitor) const
    {
      GougeBlock3D::visitParticles(visitor);
      for(vector<SBG_ptr>::iterator it = m_grainParticleGen.begin();
	  it != m_grainParticleGen.end();
	  it++)
	{
	  SphereBlockGenerator::ParticleIterator particleIt = (*it)->getParticleIterator();
          while (particleIt.hasNext()) {
            particleIt.next()->visit(visitor);
          }
	}
    }
  } // namespace lsm
} // namespace esys 
