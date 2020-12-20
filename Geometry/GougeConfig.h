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


#ifndef ESYS_LSMGOUGECONFIG_H
#define ESYS_LSMGOUGECONFIG_H

#include "Foundation/BoundingBox.h"
#include "Geometry/CircularNeighbourTable.h"
#include "Geometry/Plane3D.h"
#include "Geometry/CubicBoxPacker.h"
#include "Geometry/RandomBoxPacker.h"
#include "Geometry/PackerGenerators.h"
#include "Geometry/DistConnections.h"

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

    class ParticleRndPackPrms
    {
    public:
      inline ParticleRndPackPrms();

      inline ParticleRndPackPrms(double size, double minRadius, double maxRadius);

      inline ~ParticleRndPackPrms();

      inline double getSize() const;

      inline double getMinParticleRadius() const;

      inline double getMaxParticleRadius() const;

    private:
      double m_size;
      double m_minParticleRadius;
      double m_maxParticleRadius;
    };

    template <typename TmplParticleGrainGen>
    class GrainRndPackPrms : public ParticleRndPackPrms
    {
    public:
      typedef TmplParticleGrainGen ParticleGrainGen;
      typedef ParticleRndPackPrms  Inherited;

      GrainRndPackPrms();

      GrainRndPackPrms(
        double size,
        ParticleGrainGen &particleGrainGen,
        int connectionTag=0
      );

      double getMinGrainRadius();

      double getMaxGrainRadius();

      ParticleGrainGen &getParticleGrainGen() const;

      int getConnectionTag() const;
    private:
      ParticleGrainGen *m_pParticleGrainGen;
      int              m_connectionTag;
    };

    typedef std::vector<bool> BoolVector;
    typedef std::vector<BoundingBox> BoundingBoxVector;
    
    class PackingInfo
    {
    public:
      inline PackingInfo(
        const BoundingBox &bBox,
        const BoolVector  &periodicDimensions,
        Orientation       orientation,
        double            minRadius,
        double            maxRadius
      );
      
      inline bool is3d() const;
      
      inline void initialiseFitPlaneVector();
      
      inline const BoundingBox &getBBox() const;
      
      inline const PlaneVector &getFitPlaneVector() const;
      
      inline double getMinParticleRadius() const;
      
      inline double getMaxParticleRadius() const;
      
      inline const BoolVector &getPeriodicDimensions() const;
    private:
      BoundingBox m_bBox;
      BoolVector  m_periodicDimensions;
      Orientation m_orientation;
      double      m_minRadius;
      double      m_maxRadius;
      PlaneVector m_fitPlaneVector;
    };

    template <typename TmplParticleGrainGen>
    class GougePackingInfo : public PackingInfo
    {
    public:
      typedef TmplParticleGrainGen ParticleGrainGen;
      typedef PackingInfo          Inherited;

      GougePackingInfo(
        const BoundingBox &bBox,
        const BoolVector  &periodicDimensions,
        Orientation       orientation,
        ParticleGrainGen  &particleGrainGen
      );

      double getMinGrainRadius() const;

      double getMaxGrainRadius() const;

      ParticleGrainGen &getParticleGrainGen() const;

    private:
      ParticleGrainGen *m_pParticleGrainGen;
    };

    typedef std::vector<PackingInfo>      PackingInfoVector;

    template <typename TmplParticleGrainGen>
    class GougeConfigPrms
    {
    public:
      typedef TmplParticleGrainGen                 ParticleGrainGen;
      typedef GrainRndPackPrms<ParticleGrainGen>   GrainRPackPrms;
      typedef typename GrainRPackPrms::Inherited ParticleRndPackPrms;
      typedef GougePackingInfo<ParticleGrainGen>   GougePackInfo;
      typedef std::vector<GougePackInfo>        GougePackingInfoVector;
      
      GougeConfigPrms();
      /**
       *
       */
      GougeConfigPrms(
        const BoundingBox         &bBox,
        double                    padRadius,
        Orientation               orientation,
        const ParticleRndPackPrms &faultRegionPrms,
        const GrainRPackPrms    &gougeRegionPrms,
        const BoolVector          &peridicDimensions=BoolVector(3, false),
        int                       maxInsertionFailures=100,
        double                    tolerance = DBL_EPSILON*128,
        double                    connectionTolerance = DBL_EPSILON*128*10,
        int                       blockConnectionTag = 0
      );

      ~GougeConfigPrms();

      double getTolerance() const;

      double getConnectionTolerance() const;

      const BoundingBox &getBBox() const;

      int getMaxInsertionFailures() const;

      double getRegularBlockRadius() const;

      double getFaultMinRadius() const;

      double getFaultMaxRadius() const;

      double getGougeMinRadius() const;

      double getGougeMaxRadius() const;

      int getGougeConnectionTag() const;

      int getBlockConnectionTag() const;

      const BoolVector &getPeriodicDimensions() const;

      BoundingBoxVector getRegularBBoxVector() const;

      GougePackingInfoVector getGougePackingInfoVector() const;

      PackingInfoVector getFaultPackingInfoVector() const;

      BoundingBox cutFromCentre(double d1, double d2) const;

      Orientation getOrientation() const;

      int getOrientationIndex() const;

      double getOrientationSize() const;

      double getMaxRadius() const;
      
      double getMinRadius() const;
      
      bool is2d() const;

    private:
      BoundingBox         m_bBox;
      double              m_padRadius;
      Orientation         m_orientation;
      ParticleRndPackPrms m_faultPrms;
      GrainRPackPrms    m_gougePrms;
      BoolVector          m_periodicDimensions;
      int                 m_maxInsertionFailures;
      double              m_tolerance;
      double              m_connectionTolerance;
      int                 m_blockConnectionTag;
    };

    /**
     *
     */
    template <
      typename TmplGrainRandomBoxPacker,
      typename TmplParticleRandomBoxPacker,
      typename TmplConnection
    >
    class GougeConfig
    {
    public:
      typedef TmplConnection                                  Connection;
      typedef TmplGrainRandomBoxPacker                        GrainRandomPacker;
      typedef boost::shared_ptr<GrainRandomPacker>            GrainRandomPackerPtr;
      typedef typename GrainRandomPacker::ParticleGrainGen    ParticleGrainGen;
      typedef GougeConfigPrms<ParticleGrainGen>               GougeConfPrms;
      typedef typename GougeConfPrms::GougePackingInfoVector GougePackingInfoVector;
      typedef typename GrainRandomPacker::Particle            Particle;
      typedef typename GrainRandomPacker::ParticleIterator    ParticleIterator;
      typedef typename GrainRandomPacker::ParticleConstIterator ParticleConstIterator;
      typedef typename GrainRandomPacker::ParticleCollection  ParticleCollection;
      typedef typename GrainRandomPacker::Grain               Grain;
      typedef typename GrainRandomPacker::GrainIterator       GrainIterator;
      typedef typename GrainRandomPacker::GrainConstIterator  GrainConstIterator;
      typedef typename GrainRandomPacker::GrainCollection     GrainCollection;

      typedef typename GrainRandomPacker::PackerBase          APacker;
      typedef typename GrainRandomPacker::BoxPackerBase       ABoxPacker;
      
      typedef ConstRadiusGen<Particle>                        RegRadiusGen;
      typedef CubicBoxPacker<RegRadiusGen,ABoxPacker>         RegBoxPacker;
      typedef typename RegBoxPacker::ParticleGeneratorPtr     RegRadiusGenPtr;

      typedef TmplParticleRandomBoxPacker                     RndBoxPacker;
      typedef typename RndBoxPacker::ParticleGenerator        RndRadiusGen;
      typedef typename RndBoxPacker::ParticleGeneratorPtr     RndRadiusGenPtr;

      typedef typename GrainRandomPacker::NTable              NTable;
      typedef typename GrainRandomPacker::NTablePtr           NTablePtr;
      typedef boost::shared_ptr<APacker>                      GeneratorPtr;
      typedef std::vector<GeneratorPtr>                       GeneratorPtrVector;
      typedef std::vector<GrainRandomPackerPtr>               GrainRndPackerPtrVector;
      typedef typename GrainRandomPacker::ParticlePool        ParticlePool;
      typedef typename GrainRandomPacker::ParticlePoolPtr     ParticlePoolPtr;
      typedef typename GrainRandomPacker::GrainPool           GrainPool;
      typedef typename GrainRandomPacker::GrainPoolPtr        GrainPoolPtr;

      class ConnectionCmp
      {
      public:
        bool operator()(const Connection &i1, const Connection &i2) const
        {
          return
            (
              (i1.getP1Id() < i2.getP1Id())
              ||
              (
                (i1.getP1Id() == i2.getP1Id())
                &&
                (
                  (i1.getP2Id() < i2.getP2Id())
                  ||
                  (
                    (i1.getP2Id() == i2.getP2Id())
                    &&
                    (i1.getTag() < i2.getTag())
                  )
                )
              )
            );
        }

        bool operator()(const Connection *i1, const Connection *i2) const
        {
          return (*this)(*i1, *i2);
        }
      };
      typedef std::set<Connection,ConnectionCmp>  ConnectionSet;
      typedef DistConnections<Particle,Connection> ConnectionFinder;
      
      GougeConfig(const GougeConfPrms &prms);

      virtual ~GougeConfig();

      virtual void generate();

      int getNumParticles() const;

      int getNumGrains() const;

      int getNumConnections() const;

      const GrainRndPackerPtrVector &getGougeGeneratorVector() const;

      GrainRndPackerPtrVector &getGougeGeneratorVector();

      const GeneratorPtrVector &getFaultGeneratorVector() const;

      bool isGougeParticle(const Particle &particle) const;

      bool areInDifferentFaultBlocks(
        const Particle &p1,
        const Particle &p2
      ) const;

      virtual void write(std::ostream &oStream) const;

      void writeToFile(const std::string &fileName) const;

      void tagGougeParticles(int tag);

      void tagRndBlockParticles(int tag);
      
      void tagDrivingPlateParticles(
        int minDrivingTag,
        int maxDrivingTag,
        double distanceFromBBoxEdge
      );

      virtual void createConnectionSet();

      const ConnectionSet &getConnectionSet() const;

      GrainCollection getGrainCollection();

      ParticleCollection getParticleCollection();

      template <typename TmplVisitor>
      void visitParticles(TmplVisitor &visitor)
      {
        for (
          typename GeneratorPtrVector::iterator it = m_genPtrVector.begin();
          it != m_genPtrVector.end();
          it++
        )
        {
          ParticleIterator particleIt = (*it)->getParticleIterator();
          while (particleIt.hasNext()) {
            particleIt.next().visit(visitor);
          }
        }
      }

      template <typename TmplVisitor>
      void visitParticles(const TmplVisitor &visitor) const
      {
        for (
          typename GeneratorPtrVector::const_iterator it = m_genPtrVector.begin();
          it != m_genPtrVector.end();
          it++
        )
        {
          ParticleIterator particleIt = (*it)->getParticleIterator();
          while (particleIt.hasNext()) {
            particleIt.next().visit(visitor);
          }
        }
      }

      template <typename TmplVisitor>
      void visitConnections(TmplVisitor &visitor) const
      {
        const ConnectionSet &connectionSet = getConnectionSet();
        for (
          typename ConnectionSet::const_iterator it = connectionSet.begin();
          it != connectionSet.end();
          it++
        )
        {
          it->visit(visitor);
        }
      }
      
      const GougeConfPrms &getPrms() const
      {
        return m_prms;
      }

      class IdCompare
      {
      public:
        bool operator()(const Particle *p1, const Particle *p2) const
        {
          return (p1->getID() < p2->getID());
        }
      };

      class ConnectionValidator
      {
      public:
        inline ConnectionValidator(const GougeConfig &gougeBlock, double tolerance)
          : m_pGougeConfig(&gougeBlock),
            m_tolerance(tolerance)
        {
        }
  
        inline bool isValid(const Particle &p1, const Particle &p2) const
        {
          return
            (
              (p1.getID() < p2.getID())
              &&
              ((p1.getPos() - p2.getPos()).norm() < (m_tolerance + (p1.getRad() + p2.getRad())))
              &&
              ((!m_pGougeConfig->isGougeParticle(p1)) && (!m_pGougeConfig->isGougeParticle(p2)))
              &&
              ((!m_pGougeConfig->areInDifferentFaultBlocks(p1, p2)))
            );
        }
  
      private:
        const GougeConfig *m_pGougeConfig;
        double             m_tolerance;
      };
  
      class GeoParticleWriter
      {
      public:
        GeoParticleWriter(std::ostream &oStream, int precision)
          : m_pOStream(&oStream),
            m_precision(precision)
        {
        }
        
        void visitParticle(const Particle &particle) const
        {
          (*m_pOStream)
            << std::setprecision(m_precision)
            << particle.getPos() << " "
            << particle.getRad() << " "
            << particle.getID() << " "
            << particle.getTag() << "\n";
        }
  
      private:
        std::ostream *m_pOStream;
        int          m_precision;
      };
  
      class GeoConnectionWriter
      {
      public:
        GeoConnectionWriter(std::ostream &oStream)
          : m_pOStream(&oStream)
        {
        }
  
        void visitBasicInteraction(const BasicInteraction &connection)
        {
          (*m_pOStream)
            << connection.first()  << " "
            << connection.second() << " "
            << 0 << "\n";
        }
  
      private:
        std::ostream *m_pOStream;
        int          m_precision;
      };

    protected:
      NTablePtr               m_nTablePtr;
      GougeConfPrms         m_prms;
      ConnectionSet           m_connectionSet;
      GrainRndPackerPtrVector m_gougeGenPtrVector;
      GeneratorPtrVector      m_genPtrVector;
      ParticlePoolPtr         m_particlePoolPtr;
      GrainPoolPtr            m_grainPoolPtr;

      void createRegularBlockGenerators();
      void createFaultBlockGenerators();
      virtual void createGougeConfigGenerators();

    private:
      GeneratorPtrVector m_regularGenPtrVector;
      GeneratorPtrVector m_faultGenPtrVector;
    };
  }
}

#include "Geometry/GougeConfig.hpp"

#endif
