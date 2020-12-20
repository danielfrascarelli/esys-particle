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


#ifndef ESYS_LSMGRAINRANDOMBOXPACKER_H
#define ESYS_LSMGRAINRANDOMBOXPACKER_H

#include "Geometry/RandomBoxPacker.h"
#include "Geometry/GrainCollection.h"

#include <vector>
#include <boost/shared_ptr.hpp>

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    template <typename TmplParticleGrainGen, typename TmplGrainCollection, template <typename TPartGrainGen>  class TmplRndBoxPackerWrap>
    class GrainRandomBoxPacker : public TmplRndBoxPackerWrap<TmplParticleGrainGen>::RandomBoxPackerBase
    {
    public:
      typedef TmplGrainCollection                          GrainCollection;
      typedef boost::shared_ptr<GrainCollection>           GrainCollectionPtr;
      typedef typename GrainCollection::GrainIterator      GrainIterator;
      typedef typename GrainCollection::GrainConstIterator GrainConstIterator;
      typedef typename GrainCollection::Grain              Grain;
      typedef typename Grain::Id                           GrainId;
      typedef typename GrainCollection::GrainPool          GrainPool;
      typedef typename GrainCollection::GrainPoolPtr       GrainPoolPtr;
      typedef
        typename TmplRndBoxPackerWrap<TmplParticleGrainGen>::RandomBoxPackerBase
        Inherited;
      typedef Inherited                                    RandomBoxPackerBase;
      typedef typename Inherited::ParticleGenerator        ParticleGrainGen;
      typedef typename Inherited::ParticleGeneratorPtr     ParticleGrainGenPtr;
      typedef typename Inherited::Particle                 Particle;
      typedef typename Inherited::NTable                   NTable;
      typedef typename Inherited::NTablePtr                NTablePtr;
      typedef typename Inherited::ParticleVector           ParticleVector;
      typedef typename Inherited::ParticlePool             ParticlePool;
      typedef typename Inherited::ParticlePoolPtr          ParticlePoolPtr;
      typedef typename Inherited::PlaneVector              PlaneVector;
      typedef typename Inherited::StuffedParticleIterator  StuffedParticleIterator;

      GrainRandomBoxPacker(
        ParticleGrainGenPtr particleGrainGenPtr,
        ParticlePoolPtr     particlePoolPtr,
        NTablePtr           nTablePtr,
        const BoundingBox   &bBox,
        const BoolVector    &periodicDimensions,
        double              tolerance,
        double              cubicPackRadius,
        int                 maxInsertionFailures
      );

      GrainRandomBoxPacker(
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
      );

      virtual ~GrainRandomBoxPacker();

      ParticleGrainGen &getParticleGrainGen();

      const ParticleGrainGen &getParticleGrainGen() const;

      void setParticleGrainGen(ParticleGrainGen &particleGrainGen);

      void setParticleGrainGen(ParticleGrainGenPtr particleGrainGenPtr);

      GrainId getNextGrainId() const;
      
      Grain &constructGrain();

      Grain &createAndInsertGrain(const Grain &grain);

      void generateRandomFillGrains();

      void generateCubicPackingGrains();

      int getNumGrains() const;

      GrainConstIterator getGrainIterator() const;

      GrainIterator getGrainIterator();

      const GrainCollection &getGrainCollection() const;

      GrainCollection &getGrainCollection();

      virtual void generate();

    private:
      GrainCollectionPtr m_grainCollectionPtr;
    };
  }
}

#include "Geometry/GrainRandomBoxPacker.hpp"

#endif
