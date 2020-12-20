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


#include "Foundation/console.h"
#include "Geometry/GrainRandomBoxPacker.h"
#include "Geometry/SphereFitter.h"

#include <algorithm>
#include <stdexcept>
#include <float.h>

namespace esys
{
  namespace lsm
  {
    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::GrainRandomBoxPacker(
      ParticleGrainGenPtr particleGrainGenPtr,
      ParticlePoolPtr     particlePoolPtr,
      NTablePtr           nTablePtr,
      const BoundingBox   &bBox,
      const BoolVector    &periodicDimensions,
      double              tolerance,
      double              cubicPackRadius,
      int                 maxInsertionFailures,
      const PlaneVector   &fitPlaneVector,
      GrainPoolPtr        grainPoolPtr
    )
     :  Inherited(
          particleGrainGenPtr,
          particlePoolPtr,
          nTablePtr,
          bBox,
          periodicDimensions,
          tolerance,
          cubicPackRadius,
          maxInsertionFailures,
          fitPlaneVector
        ),
        m_grainCollectionPtr(new GrainCollection(particlePoolPtr, grainPoolPtr))
    {
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::GrainRandomBoxPacker(
      ParticleGrainGenPtr particleGrainGenPtr,
      ParticlePoolPtr     particlePoolPtr,
      NTablePtr           nTablePtr,
      const BoundingBox   &bBox,
      const BoolVector    &periodicDimensions,
      double              tolerance,
      double              cubicPackRadius,
      int                 maxInsertionFailures
    )
     :  Inherited(
          particleGrainGenPtr,
          particlePoolPtr,
          nTablePtr,
          bBox,
          periodicDimensions,
          tolerance,
          cubicPackRadius,
          maxInsertionFailures
        ),
        m_grainCollectionPtr(new GrainCollection(particlePoolPtr))
    {
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::~GrainRandomBoxPacker()
    {
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    const typename GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::GrainCollection &
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::getGrainCollection() const
    {
      return *(m_grainCollectionPtr.get());
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    typename GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::GrainCollection &
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::getGrainCollection()
    {
      return *(m_grainCollectionPtr.get());
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    typename GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::GrainIterator
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::getGrainIterator()
    {
      return getGrainCollection().getGrainIterator();
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    typename GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::GrainConstIterator
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::getGrainIterator() const
    {
      return getGrainCollection().getGrainIterator();
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    int
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::getNumGrains() const
    {
      return getGrainCollection().getNumGrains();
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    typename GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::ParticleGrainGen &
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::getParticleGrainGen()
    {
      return Inherited::getParticleGenerator();
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    const typename GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::ParticleGrainGen &
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::getParticleGrainGen() const
    {
      return Inherited::getParticleGenerator();
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    void
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::setParticleGrainGen(
      ParticleGrainGen &particleGrainGen
    )
    {
      this->setParticleGenerator(particleGrainGen);
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    void
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::setParticleGrainGen(
      ParticleGrainGenPtr particleGrainGenPtr
    )
    {
      setParticleGenerator(particleGrainGenPtr);
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    typename GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::GrainId
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::getNextGrainId() const
    {
      return m_grainCollectionPtr->getNumGrains();
    }
    
    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    typename GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::Grain &
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::constructGrain()
    {
      return m_grainCollectionPtr->createGrain();
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    typename GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::Grain &
    GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::createAndInsertGrain(
      const Grain &grain
    )
    {
      Grain &g = constructGrain();
      g.setId(getNextGrainId());
      typename Grain::ParticleConstIterator it = grain.getParticleIterator();
      while (it.hasNext())
      {
        g.insertRef(this->createAndInsertParticle(it.next()));
      }
      return g;
    }
    
    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    void GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::generateCubicPackingGrains()
    {
      GridIterator pointIt = GridIterator(this->getBBox(), this->getCubicPackingRadius());
      while (pointIt.hasNext()) {
        const Particle candidate =
          this->getCandidateParticle(pointIt.next(), this->getCubicPackingRadius());
        if (this->particleFitsInBBoxWithNeighbours(candidate)) {
          createAndInsertGrain(getParticleGrainGen().getGrain(candidate));
        }
      }
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    void GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::generateRandomFillGrains()
    {
      StuffedParticleIterator it =
        StuffedParticleIterator(
          *this,
          this->getMaxInsertionFailures(),
          this->getFitPlaneVector()
        );
      while (it.hasNext())
      {
        createAndInsertGrain(getParticleGrainGen().getGrain(it.next()));
      }
      it.logInfo();
    }

    template <typename TGrainGen, typename TGrainCol, template <typename TTGrainGen> class TRndPackWrap>
    void GrainRandomBoxPacker<TGrainGen,TGrainCol,TRndPackWrap>::generate()
    {
      generateCubicPackingGrains();
      generateRandomFillGrains();
    }
  }
}
