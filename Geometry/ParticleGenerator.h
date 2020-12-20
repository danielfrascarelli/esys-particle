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


#ifndef ESYS_LSMPARTICLEGENERATOR_H
#define ESYS_LSMPARTICLEGENERATOR_H

#include "Geometry/SimpleParticle.h"
#include "Geometry/CircularNeighbourTable.h"

#include <boost/pool/object_pool.hpp>

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    class ParticleGenerator
    {
    public:
      typedef CircularNeighbourTable<SimpleParticle> NTable;
      typedef boost::object_pool<SimpleParticle>     ParticlePool;
      
      ParticleGenerator(NTable &nTable, ParticlePool &particlePool);

      virtual ~ParticleGenerator();

      virtual void generate() = 0;
      
    protected:
      ParticleGenerator();
      
      NTable &getNTable();
      const NTable &getNTable() const;
      
      ParticlePool &getParticlePool();
      const ParticlePool &getParticlePool() const;
    private:
      NTable       *m_pNTable;
      ParticlePool *m_pParticlePool;
    };
  }
}

#endif
