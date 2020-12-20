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


#include "Geometry/ParticleGenerator.h"
#include "Geometry/CircularNeighbourTable.h"

namespace esys
{
  namespace lsm
  {
    ParticleGenerator::ParticleGenerator()
      : m_pNTable(NULL),
        m_pParticlePool(NULL)
    {
    }

    ParticleGenerator::ParticleGenerator(
      NTable       &nTable,
      ParticlePool &particlePool
    )
      : m_pNTable(&nTable),
        m_pParticlePool(&particlePool)
    {
    }

    ParticleGenerator::~ParticleGenerator()
    {
    }
    
    ParticleGenerator::NTable &ParticleGenerator::getNTable()
    {
      return *(m_pNTable);
    }

    const ParticleGenerator::NTable &ParticleGenerator::getNTable() const
    {
      return *(m_pNTable);
    }

    ParticleGenerator::ParticlePool &ParticleGenerator::getParticlePool()
    {
      return *m_pParticlePool;
    }
    
    const ParticleGenerator::ParticlePool &ParticleGenerator::getParticlePool() const
    {
      return *m_pParticlePool;
    }
  }
}
