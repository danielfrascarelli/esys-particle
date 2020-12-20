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


#ifndef ESYS_LSMGOUGEBLOCK3D_H
#define ESYS_LSMGOUGEBLOCK3D_H

#include "Foundation/BoundingBox.h"
#include "Geometry/CircularNeighbourTable.h"
#include "Geometry/BlockGenerator.h"
#include "Geometry/Plane3D.h"
#include "Geometry/BasicInteraction.h"

#include <boost/shared_ptr.hpp>

#include <vector>
#include <float.h>

namespace esys
{
  namespace lsm
  {
    typedef std::vector<Plane3D> PlaneVector;
    enum Orientation
    {
      XY,
      XZ,
      YZ
    };

    class ParticleBlockPrms
    {
    public:
      ParticleBlockPrms();

      ParticleBlockPrms(double size, double minRadius, double maxRadius);

      ~ParticleBlockPrms();

      double m_size;
      double m_minRadius;
      double m_maxRadius;
    };

    typedef std::vector<bool> BoolVector;
    typedef std::vector<BoundingBox> BoundingBoxVector;
    
    class PackingInfo
    {
    public:
      PackingInfo(
        const BoundingBox &bBox,
        const BoolVector  &periodicDimensions,
        Orientation       orientation,
        double            minRadius,
        double            maxRadius
      );
      
      bool is3d() const;
      
      void initialiseFitPlaneVector();
      
      const BoundingBox &getBBox() const;
      
      const PlaneVector &getFitPlaneVector() const;
      
      double getMinRadius() const;
      
      double getMaxRadius() const;
      
      const BoolVector &getPeriodicDimensions() const;
    private:
      BoundingBox m_bBox;
      BoolVector  m_periodicDimensions;
      Orientation m_orientation;
      double      m_minRadius;
      double      m_maxRadius;
      PlaneVector m_fitPlaneVector;
    };

    typedef std::vector<PackingInfo> PackingInfoVector;
    class GougeBlockPrms
    {
    public:
      GougeBlockPrms();
      /**
       *
       */
      GougeBlockPrms(
        const BoundingBox       &bBox,
        double                  padRadius,
        Orientation             orientation,
        const ParticleBlockPrms &faultRegionPrms,
        const ParticleBlockPrms &gougeRegionPrms,
        const BoolVector        &peridicDimensions=BoolVector(3, false),
        int                     maxInsertionFailures=100,
        double                  tolerance = DBL_EPSILON*128,
        double                  connectionTolerance = DBL_EPSILON*128*10
      );

      ~GougeBlockPrms();

      double getTolerance() const;

      double getConnectionTolerance() const;

      const BoundingBox &getBBox() const;

      int getMaxInsertionFailures() const;

      double getRegularBlockRadius() const;

      double getFaultMinRadius() const;

      double getFaultMaxRadius() const;

      double getGougeMinRadius() const;

      double getGougeMaxRadius() const;

      const BoolVector &getPeriodicDimensions() const;

      BoundingBoxVector getRegularBBoxVector() const;

      PackingInfoVector getGougePackingInfoVector() const;

      PackingInfoVector getFaultPackingInfoVector() const;

      BoundingBox cutFromCentre(double d1, double d2) const;

      Orientation getOrientation() const;

      int getOrientationIndex() const;

      double getOrientationSize() const;

      double getMaxRadius() const;
      
      double getMinRadius() const;
      
      bool is2d() const;

    private:
      BoundingBox       m_bBox;
      double            m_padRadius;
      Orientation       m_orientation;
      ParticleBlockPrms m_faultPrms;
      ParticleBlockPrms m_gougePrms;
      BoolVector        m_periodicDimensions;
      int               m_maxInsertionFailures;
      double            m_tolerance;
      double            m_connectionTolerance;
    };

    /*!
      \class GougeBlock3D
      \brief Block consisting of regular padding, random layer and gouge
     */
    class GougeBlock3D
    {
    public:
      typedef SimpleParticle Particle;
      GougeBlock3D(const GougeBlockPrms &prms);

      virtual ~GougeBlock3D();

      virtual void generate();

      int getNumParticles() const;

      typedef CircularNeighbourTable<SimpleParticle> NTable;
      typedef boost::shared_ptr<NTable>              NTablePtr;
      typedef boost::shared_ptr<BlockGenerator>      GeneratorPtr;
      typedef std::vector<GeneratorPtr>              GeneratorPtrVector;
      typedef NTable::ParticlePool                   ParticlePool;
      typedef NTable::ParticlePoolPtr                ParticlePoolPtr;

      const GeneratorPtrVector &getGougeGeneratorVector() const;

      const GeneratorPtrVector &getFaultGeneratorVector() const;

      bool isGougeParticle(const SimpleParticle &particle) const;

      bool areInDifferentFaultBlocks(
        const SimpleParticle &p1,
        const SimpleParticle &p2
      ) const;

      virtual void write(std::ostream &oStream) const;

      void writeToFile(const std::string &fileName) const;

      void tagGougeParticles(int tag);

      void tagFaultParticles(int tag);
      
      void tagDrivingPlateParticles(
        int minDrivingTag,
        int maxDrivingTag,
        double distanceFromBBoxEdge
      );

      typedef std::set<BasicInteraction,BILess> InteractionSet;
      
      virtual void createInteractionSet();

      const InteractionSet &getInteractionSet() const;

      template <typename TmplVisitor>
      void visitParticles(TmplVisitor &visitor)
      {
        for (
          GeneratorPtrVector::iterator it = m_genPtrVector.begin();
          it != m_genPtrVector.end();
          it++
        )
        {
          BlockGenerator::ParticleIterator particleIt = (*it)->getParticleIterator();
          while (particleIt.hasNext()) {
            particleIt.next()->visit(visitor);
          }
        }
      }

      template <typename TmplVisitor>
      void visitParticles(const TmplVisitor &visitor) const
      {
        for (
          GeneratorPtrVector::const_iterator it = m_genPtrVector.begin();
          it != m_genPtrVector.end();
          it++
        )
        {
          BlockGenerator::ParticleIterator particleIt = (*it)->getParticleIterator();
          while (particleIt.hasNext()) {
            particleIt.next()->visit(visitor);
          }
        }
      }

      template <typename TmplVisitor>
      void visitInteractions(TmplVisitor &visitor) const
      {
        const InteractionSet &interactionSet = getInteractionSet();
        for (
          InteractionSet::const_iterator it = interactionSet.begin();
          it != interactionSet.end();
          it++
        )
        {
          it->visit(visitor);
        }
      }
      
      const GougeBlockPrms &getPrms() const
      {
        return m_prms;
      }

    protected:
      NTablePtr          m_nTablePtr;
      GougeBlockPrms     m_prms;
      InteractionSet     m_interactionSet;
      GeneratorPtrVector m_gougeGenPtrVector;
      GeneratorPtrVector m_genPtrVector;
      ParticlePoolPtr    m_particlePoolPtr;

      void createRegularBlockGenerators();
      void createFaultBlockGenerators();
      virtual void createGougeBlockGenerators();

    private:
      GeneratorPtrVector m_regularGenPtrVector;
      GeneratorPtrVector m_faultGenPtrVector;
    };
  }
}

#endif
