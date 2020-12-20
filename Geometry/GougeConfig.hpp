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
#include "Geometry/GeometryInfo.h"
#include "Geometry/CircularNeighbourTable.h"

#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace esys
{
  namespace lsm
  {
    //==========================================================================
    PackingInfo::PackingInfo(
      const BoundingBox &bBox,
      const BoolVector  &periodicDimensions,
      Orientation       orientation,
      double            minRadius,
      double            maxRadius
    ) : m_bBox(bBox),
        m_periodicDimensions(periodicDimensions),
        m_orientation(orientation),
        m_minRadius(minRadius),
        m_maxRadius(maxRadius)
    {
      initialiseFitPlaneVector();
    }
    
    bool PackingInfo::is3d() const
    {
      return (m_bBox.getSizes().Z() > 0.0);
    }
    
    void PackingInfo::initialiseFitPlaneVector()
    {
      m_fitPlaneVector.clear();
      if ((m_orientation != XZ) && (!getPeriodicDimensions()[1])) {
        m_fitPlaneVector.push_back(
          Plane3D(Vec3(0,  1, 0), getBBox().getMinPt())
        );
        m_fitPlaneVector.push_back(
          Plane3D(Vec3(0, -1, 0), getBBox().getMaxPt())
        );
      }
      if ((m_orientation != YZ) && (!getPeriodicDimensions()[0])) {
        m_fitPlaneVector.push_back(
          Plane3D(Vec3( 1, 0, 0), getBBox().getMinPt())
        );
        m_fitPlaneVector.push_back(
          Plane3D(Vec3(-1, 0, 0), getBBox().getMaxPt())
        );
      }
      if (
          is3d()
          &&
          (m_orientation != XY)
          &&
          (!getPeriodicDimensions()[2])
      ) {
        m_fitPlaneVector.push_back(
          Plane3D(Vec3(0, 0,  1), getBBox().getMinPt())
        );
        m_fitPlaneVector.push_back(
          Plane3D(Vec3(0, 0, -1), getBBox().getMaxPt())
        );
      }
    }

    const BoundingBox &PackingInfo::getBBox() const
    {
      return m_bBox;
    }

    const PlaneVector &PackingInfo::getFitPlaneVector() const
    {
      return m_fitPlaneVector;
    }

    double PackingInfo::getMinParticleRadius() const
    {
      return m_minRadius;
    }

    double PackingInfo::getMaxParticleRadius() const
    {
      return m_maxRadius;
    }

    const BoolVector &PackingInfo::getPeriodicDimensions() const
    {
      return m_periodicDimensions;
    }

    //==========================================================================
    
    template <typename TGrainGen>
    GougePackingInfo<TGrainGen>::GougePackingInfo(
      const BoundingBox &bBox,
      const BoolVector  &periodicDimensions,
      Orientation       orientation,
      ParticleGrainGen  &particleGrainGen
    ) : Inherited(
          bBox,
          periodicDimensions,
          orientation,
          particleGrainGen.getMinParticleRadius(),
          particleGrainGen.getMaxParticleRadius()
        ),
        m_pParticleGrainGen(&particleGrainGen)
        
    {
    }

    template <typename TGrainGen>
    typename GougePackingInfo<TGrainGen>::ParticleGrainGen &
    GougePackingInfo<TGrainGen>::getParticleGrainGen() const
    {
      return *m_pParticleGrainGen;
    }

#if 0
    template <typename TGrainGen>
    const typename GougePackingInfo<TGrainGen>::ParticleGrainGen &
    GougePackingInfo<TGrainGen>::getParticleGrainGen() const
    {
      return *m_pParticleGrainGen;
    }
#endif

    template <typename TGrainGen>
    double GougePackingInfo<TGrainGen>::getMinGrainRadius() const
    {
      return getParticleGrainGen().getMinGrainRadius();
    }

    template <typename TGrainGen>
    double GougePackingInfo<TGrainGen>::getMaxGrainRadius() const
    {
      return getParticleGrainGen().getMaxGrainRadius();
    }

    //==========================================================================

    ParticleRndPackPrms::ParticleRndPackPrms()
      : m_size(0.0),
        m_minParticleRadius(0.0),
        m_maxParticleRadius(0.0)
    {
    }

    ParticleRndPackPrms::ParticleRndPackPrms(
      double size,
      double minRadius,
      double maxRadius
    ) : m_size(size),
        m_minParticleRadius(minRadius),
        m_maxParticleRadius(maxRadius)
    {
    }

    ParticleRndPackPrms::~ParticleRndPackPrms()
    {
    }

    double ParticleRndPackPrms::getSize() const
    {
      return m_size;
    }

    double ParticleRndPackPrms::getMinParticleRadius() const
    {
      return m_minParticleRadius;
    }

    double ParticleRndPackPrms::getMaxParticleRadius() const
    {
      return m_maxParticleRadius;
    }

    //==========================================================================
    
    template <typename TPGrainGen>
    GrainRndPackPrms<TPGrainGen>::GrainRndPackPrms()
      :
        Inherited(),
        m_pParticleGrainGen(NULL)
    {
    }

    template <typename TPGrainGen>
    GrainRndPackPrms<TPGrainGen>::GrainRndPackPrms(
      double size,
      ParticleGrainGen &particleGrainGen,
      int connectionTag
    ) : Inherited(
          size,
          particleGrainGen.getMinParticleRadius(),
          particleGrainGen.getMaxParticleRadius()
        ),
        m_pParticleGrainGen(&particleGrainGen),
        m_connectionTag(connectionTag)
    {
    }

    template <typename TPGrainGen>
    typename GrainRndPackPrms<TPGrainGen>::ParticleGrainGen &
    GrainRndPackPrms<TPGrainGen>::getParticleGrainGen() const
    {
      return *m_pParticleGrainGen;
    }

    template <typename TPGrainGen>
    int
    GrainRndPackPrms<TPGrainGen>::getConnectionTag() const
    {
      return m_connectionTag;
    }

#if 0
    template <typename TPGrainGen>
    const typename GrainRndPackPrms<TPGrainGen>::ParticleGrainGen &
    GrainRndPackPrms<TPGrainGen>::getParticleGrainGen() const
    {
      return *m_pParticleGrainGen;
    }
#endif

    template <typename TPGrainGen>
    double GrainRndPackPrms<TPGrainGen>::getMinGrainRadius()
    {
      return getParticleGrainGen().getMinGrainRadius();
    }

    template <typename TPGrainGen>
    double GrainRndPackPrms<TPGrainGen>::getMaxGrainRadius()
    {
      return getParticleGrainGen().getMaxGrainRadius();
    }

    //==========================================================================

    template <typename TPGrainGen>
    GougeConfigPrms<TPGrainGen>::GougeConfigPrms()
      : m_bBox(Vec3::ZERO, Vec3::ZERO),
        m_padRadius(0.0),
        m_orientation(XZ),
        m_faultPrms(),
        m_gougePrms(),
        m_periodicDimensions(3, false),
        m_maxInsertionFailures(50),
        m_tolerance(DBL_EPSILON*128),
        m_connectionTolerance(DBL_EPSILON*128*10),
        m_blockConnectionTag(0)
    {
    }

    template <typename TPGrainGen>
    GougeConfigPrms<TPGrainGen>::GougeConfigPrms(
      const BoundingBox         &bBox,
      double                    padRadius,
      Orientation               orientation,
      const ParticleRndPackPrms &faultRegionPrms,
      const GrainRPackPrms    &gougeRegionPrms,
      const BoolVector          &periodicDimensions,
      int                       maxInsertionFailures,
      double                    tolerance ,
      double                    connectionTolerance,
      int                       blockConnectionTag
    ) : m_bBox(bBox),
        m_padRadius(padRadius),
        m_orientation(orientation),
        m_faultPrms(faultRegionPrms),
        m_gougePrms(gougeRegionPrms),
        m_periodicDimensions(periodicDimensions),
        m_maxInsertionFailures(maxInsertionFailures),
        m_tolerance(tolerance),
        m_connectionTolerance(connectionTolerance),
        m_blockConnectionTag(blockConnectionTag)
    {
      m_bBox = GridIterator(bBox, getMaxRadius()).getSphereBBox();
    }

    template <typename TPGrainGen>
    GougeConfigPrms<TPGrainGen>::~GougeConfigPrms()
    {
    }

    template <typename TPGrainGen>
    const BoundingBox &GougeConfigPrms<TPGrainGen>::getBBox() const
    {
      return m_bBox;
    }

    template <typename TPGrainGen>
    int GougeConfigPrms<TPGrainGen>::getGougeConnectionTag() const
    {
      return m_gougePrms.getConnectionTag();
    }

    template <typename TPGrainGen>
    int GougeConfigPrms<TPGrainGen>::getBlockConnectionTag() const
    {
      return m_blockConnectionTag;
    }

    template <typename TPGrainGen>
    int GougeConfigPrms<TPGrainGen>::getMaxInsertionFailures() const
    {
      return m_maxInsertionFailures;
    }

    template <typename TPGrainGen>
    const BoolVector &GougeConfigPrms<TPGrainGen>::getPeriodicDimensions() const
    {
      return m_periodicDimensions;
    }

    template <typename TPGrainGen>
    Orientation GougeConfigPrms<TPGrainGen>::getOrientation() const
    {
      return m_orientation;
    }

    template <typename TPGrainGen>
    double GougeConfigPrms<TPGrainGen>::getTolerance() const
    {
      return m_tolerance;
    }

    template <typename TPGrainGen>
    double GougeConfigPrms<TPGrainGen>::getConnectionTolerance() const
    {
      return m_connectionTolerance;
    }

    template <typename TPGrainGen>
    int GougeConfigPrms<TPGrainGen>::getOrientationIndex() const
    {
      int idx = 0;
      switch (m_orientation) {
        case XZ:
        {
          idx = 1;
          break;
        }
        case YZ:
        {
          idx = 0;
          break;
        }
        case XY:
        {
          idx = 2;
          break;
        }
        default:
        {
          std::stringstream msg;
          msg << "Invalid orientation: " << m_orientation;
          throw std::runtime_error(msg.str());
        }
      }
      return idx;
    }

    template <typename TPGrainGen>
    BoundingBox GougeConfigPrms<TPGrainGen>::cutFromCentre(double d1, double d2) const
    {
      const int idx = getOrientationIndex();
      const BoundingBox bBox = getBBox(); 

      const double cmp1 = d1 + (bBox.getMaxPt()[idx] + bBox.getMinPt()[idx])/2.0;
      const double cmp2 = d2 + (bBox.getMaxPt()[idx] + bBox.getMinPt()[idx])/2.0;
      Vec3 minPt = bBox.getMinPt();
      Vec3 maxPt = bBox.getMaxPt();

      minPt[idx] = std::min(cmp1, cmp2);
      maxPt[idx] = std::max(cmp1, cmp2);

      Vec3 tmpPt = maxPt;
      tmpPt[idx] = minPt[idx];

      if ((tmpPt - bBox.getMinPt())[idx] > getTolerance()) {
        const BoundingBox minBBox = GridIterator(BoundingBox(bBox.getMinPt(), tmpPt), getMaxRadius()).getSphereBBox();
        tmpPt = minPt;
        tmpPt[idx] = minBBox.getMaxPt()[idx];
      }
      else {
        tmpPt = bBox.getMinPt();
      }
      const BoundingBox maxBBox = GridIterator(BoundingBox(bBox.getMinPt(), maxPt), getMaxRadius()).getSphereBBox();
      BoundingBox returnBBox = BoundingBox(tmpPt, maxBBox.getMaxPt());

      return returnBBox;
    }

    template <typename TPGrainGen>
    double GougeConfigPrms<TPGrainGen>::getRegularBlockRadius() const
    {
      return m_padRadius;
    }

    template <typename TPGrainGen>
    double GougeConfigPrms<TPGrainGen>::getFaultMinRadius() const
    {
      return m_faultPrms.getMinParticleRadius();
    }

    template <typename TPGrainGen>
    double GougeConfigPrms<TPGrainGen>::getFaultMaxRadius() const
    {
      return m_faultPrms.getMaxParticleRadius();
    }

    template <typename TPGrainGen>
    double GougeConfigPrms<TPGrainGen>::getGougeMinRadius() const
    {
      return m_gougePrms.getMinParticleRadius();
    }

    template <typename TPGrainGen>
    double GougeConfigPrms<TPGrainGen>::getGougeMaxRadius() const
    {
      return m_gougePrms.getMaxParticleRadius();
    }

    template <typename TPGrainGen>
    double GougeConfigPrms<TPGrainGen>::getOrientationSize() const
    {
      return 
        (
          getBBox().getMaxPt()[getOrientationIndex()]
          -
          getBBox().getMinPt()[getOrientationIndex()]
        );
    }

    template <typename TPGrainGen>
    BoundingBoxVector GougeConfigPrms<TPGrainGen>::getRegularBBoxVector() const
    {
      BoundingBoxVector bBoxVector;
      if (
        (getOrientationSize() - (m_gougePrms.getSize() + 2*m_faultPrms.getSize()))
        >
        2.0*m_padRadius
      )
      {
        bBoxVector.reserve(2);
        bBoxVector.push_back(
          cutFromCentre(
            -(m_gougePrms.getSize()/2.0 + m_faultPrms.getSize()),
            -getOrientationSize()/2.0
          )
        );
        bBoxVector.push_back(
          cutFromCentre(
            m_gougePrms.getSize()/2.0 + m_faultPrms.getSize(),
            getOrientationSize()/2.0
          )
        );
      }
      return bBoxVector;
    }

    template <typename TPGrainGen>
    typename GougeConfigPrms<TPGrainGen>::GougePackingInfoVector
    GougeConfigPrms<TPGrainGen>::getGougePackingInfoVector() const
    {
      GougePackingInfoVector infoVec;
      if (m_gougePrms.getSize() > 0.0) {
        Vec3 overlap = Vec3::ZERO;
        overlap[getOrientationIndex()] = m_faultPrms.getMaxParticleRadius();
        BoundingBox bBox =
          cutFromCentre(
             m_gougePrms.getSize()/2.0,
            -m_gougePrms.getSize()/2.0
          );
        
        infoVec.push_back(
          GougePackInfo(
            BoundingBox(bBox.getMinPt() - overlap, bBox.getMaxPt() + overlap),
            getPeriodicDimensions(),
            getOrientation(),
            m_gougePrms.getParticleGrainGen()
          )
        );
      }
      return infoVec;
    }

    template <typename TPGrainGen>
    PackingInfoVector GougeConfigPrms<TPGrainGen>::getFaultPackingInfoVector() const
    {
      PackingInfoVector infoVec;
      if (m_faultPrms.getSize() > 0.0)
      {
        if (
            (getOrientationSize() - (m_gougePrms.getSize() + 2.0*m_faultPrms.getSize()))
            >
            0.0
        )
        {
          infoVec.reserve(2);
          const double roughnessSize = m_faultPrms.getSize();
          Vec3 overlap = Vec3::ZERO;
          overlap[getOrientationIndex()] = m_padRadius;
          const BoundingBox bBox1 =
            cutFromCentre(
              -m_gougePrms.getSize()/2.0,
              -(m_gougePrms.getSize()/2.0 + roughnessSize)
            );
          const BoundingBox bBox2 =
            cutFromCentre(
              m_gougePrms.getSize()/2.0,
              m_gougePrms.getSize()/2.0 + roughnessSize
            );
  
          infoVec.push_back(
            PackingInfo(
              BoundingBox(bBox1.getMinPt() - overlap, bBox1.getMaxPt()),
              getPeriodicDimensions(),
              getOrientation(),
              getFaultMinRadius(),
              getFaultMaxRadius()
            )
          );
  
          infoVec.push_back(
            PackingInfo(
              BoundingBox(bBox2.getMinPt(), bBox2.getMaxPt() + overlap),
              getPeriodicDimensions(),
              getOrientation(),
              getFaultMinRadius(),
              getFaultMaxRadius()
            )
          );
        }
        else
        {
          std::stringstream msg;
          msg
            << "Roughness size plus gouge size is greater than block size: "
            << "2*" << m_faultPrms.getSize() << " + " << m_gougePrms.getSize()
            << " > " << getOrientationSize();
          throw std::runtime_error(msg.str().c_str());
        }
      }
      return infoVec;
    }

    template <typename TPGrainGen>
    double GougeConfigPrms<TPGrainGen>::getMaxRadius() const
    {
      return 
        std::max(
          m_padRadius,
          std::max(m_faultPrms.getMaxParticleRadius(), m_gougePrms.getMaxParticleRadius())
        );
    }

    template <typename TPGrainGen>
    double GougeConfigPrms<TPGrainGen>::getMinRadius() const
    {
      return 
        std::min(
          m_padRadius,
          std::min(m_faultPrms.getMinParticleRadius(), m_gougePrms.getMinParticleRadius())
        );
    }

    template <typename TPGrainGen>
    bool GougeConfigPrms<TPGrainGen>::is2d() const
    {
      return (getBBox().getSizes()[2] == 0.0);
    }

    //==========================================================================

    template <typename TGPckr,typename TPPckr,typename TConn>
    int GougeConfig<TGPckr,TPPckr,TConn>::getNumParticles() const
    {
      int numParticles = 0;
      for (
        typename GeneratorPtrVector::const_iterator it = m_genPtrVector.begin();
        it != m_genPtrVector.end();
        it++
      )
      {
        numParticles += (*it)->getNumParticles();
      }
      return numParticles;
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    int GougeConfig<TGPckr,TPPckr,TConn>::getNumGrains() const
    {
      int numGrains = 0;
      for (
        typename GrainRndPackerPtrVector::const_iterator packerIt =
          getGougeGeneratorVector().begin();
        packerIt != getGougeGeneratorVector().end();
        packerIt++
      )
      {
        numGrains += (*packerIt)->getNumGrains();
      }
      return numGrains;
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    int GougeConfig<TGPckr,TPPckr,TConn>::getNumConnections() const
    {
      return m_connectionSet.size();
    }
    
    template <typename TGPckr,typename TPPckr,typename TConn>
    GougeConfig<TGPckr,TPPckr,TConn>::GougeConfig(const GougeConfPrms &prms)
      : m_nTablePtr(),
        m_prms(prms),
        m_connectionSet(),
        m_gougeGenPtrVector(),
        m_genPtrVector(),
        m_particlePoolPtr(new ParticlePool(8*4096)),
        m_grainPoolPtr(new GrainPool(4096)),
        m_regularGenPtrVector(),
        m_faultGenPtrVector()
    {
      /*
       * Adjust the size of the ntable bounding-box to accommodate circular
       * boundary conditions.
       */
      const BoundingBox bBox = m_prms.getBBox();
      Vec3 ntableAdjust = 
        Vec3(
          m_prms.getPeriodicDimensions()[0] ? 1 : 0,
          m_prms.getPeriodicDimensions()[1] ? 1 : 0,
          m_prms.getPeriodicDimensions()[2] ? 1 : 0
        )*m_prms.getMaxRadius();
      if (m_prms.getBBox().getSizes().Z() >= 4*m_prms.getMaxRadius()) {
        ntableAdjust += Vec3(m_prms.getMaxRadius(), 0, 0);
      }
      const BoundingBox nTableBBox(bBox.getMinPt(), bBox.getMaxPt() - ntableAdjust);
      m_nTablePtr =
        NTablePtr(
          new NTable(
            nTableBBox,
            (4.0*m_prms.getMinRadius()), // grid spacing
            m_prms.getPeriodicDimensions(),
            2.1*m_prms.getMaxRadius() // width of border-region in which
                                      // particles are duplicated
                                      // for circular boundry
          )
        );
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    void GougeConfig<TGPckr,TPPckr,TConn>::createRegularBlockGenerators()
    {
      BoundingBoxVector bBoxVector = m_prms.getRegularBBoxVector();
      for (
        BoundingBoxVector::const_iterator it = bBoxVector.begin();
        it != bBoxVector.end();
        it++
      ) {
        console.Debug()
          << "GougeConfig<TGPckr,TPPckr,TConn>::createRegularBlockGenerators:"
          << "Creating RegBoxPacker in box: " << StringUtil::toString(*it)
          << "\n";

        GeneratorPtr genPtr =
          GeneratorPtr(
            new RegBoxPacker(
              RegRadiusGenPtr(new RegRadiusGen(m_prms.getRegularBlockRadius())),
              m_particlePoolPtr,
              m_nTablePtr,
              *it,
              m_prms.getPeriodicDimensions(),
              m_prms.getTolerance(),
              m_prms.getRegularBlockRadius()
            )
          );
        m_genPtrVector.push_back(genPtr);
        m_regularGenPtrVector.push_back(genPtr);
      }
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    void GougeConfig<TGPckr,TPPckr,TConn>::createFaultBlockGenerators()
    {      
      PackingInfoVector infoVec = m_prms.getFaultPackingInfoVector();
      for (
        PackingInfoVector::const_iterator it = infoVec.begin();
        it != infoVec.end();
        it++
      ) {
        console.Debug()
          << "GougeConfig<TGPckr,TPPckr,TConn>::createFaultBlockGenerators:"
          << "Creating RndBoxPacker in box: " << StringUtil::toString(it->getBBox())
          << "\n";
        GeneratorPtr genPtr = 
          GeneratorPtr(
            new RndBoxPacker(
              RndRadiusGenPtr(
                new RndRadiusGen(
                  it->getMinParticleRadius(),
                  it->getMaxParticleRadius()
                )
              ),
              m_particlePoolPtr,
              m_nTablePtr,
              it->getBBox(),
              it->getPeriodicDimensions(),
              m_prms.getTolerance(),
              it->getMaxParticleRadius(),
              m_prms.getMaxInsertionFailures(),
              it->getFitPlaneVector()
            )
          );
        m_genPtrVector.push_back(genPtr);
        m_faultGenPtrVector.push_back(genPtr);
      }
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    void GougeConfig<TGPckr,TPPckr,TConn>::createGougeConfigGenerators()
    {
      GougePackingInfoVector infoVec = m_prms.getGougePackingInfoVector();
      for (
        typename GougePackingInfoVector::const_iterator it = infoVec.begin();
        it != infoVec.end();
        it++
      ) {
        console.Debug()
          << "GougeConfig<TGPckr,TPPckr,TConn>::createGougeConfigGenerators:"
          << "Creating GrainRandomPacker in box: " << StringUtil::toString(it->getBBox())
          << "\n";
        GrainRandomPackerPtr genPtr =
          GrainRandomPackerPtr(
            new GrainRandomPacker(
              typename GrainRandomPacker::ParticleGrainGenPtr(),
              m_particlePoolPtr,
              m_nTablePtr,
              it->getBBox(),
              it->getPeriodicDimensions(),
              m_prms.getTolerance(),
              it->getParticleGrainGen().getMaxGrainRadius(),
              m_prms.getMaxInsertionFailures(),
              it->getFitPlaneVector(),
              m_grainPoolPtr
            )
          );
        genPtr->setParticleGrainGen(it->getParticleGrainGen());
        m_genPtrVector.push_back(genPtr);
        m_gougeGenPtrVector.push_back(genPtr);
      }
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    GougeConfig<TGPckr,TPPckr,TConn>::~GougeConfig()
    {
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    void GougeConfig<TGPckr,TPPckr,TConn>::generate()
    {
      // setup block generators
      createRegularBlockGenerators();
      createFaultBlockGenerators();
      createGougeConfigGenerators();
 
      // use block generators  
      console.Info() << "bbox = " << m_prms.getBBox().getMinPt() << " " << m_prms.getBBox().getMaxPt() << "\n";
      for (
        typename GeneratorPtrVector::iterator it = m_genPtrVector.begin();
        it != m_genPtrVector.end();
        it++
      )
      {
        (*it)->generate();
      }
      createConnectionSet();
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    void GougeConfig<TGPckr,TPPckr,TConn>::writeToFile(const std::string &fileName) const
    {
      std::ofstream fStream(fileName.c_str());
      write(fStream);
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    const typename GougeConfig<TGPckr,TPPckr,TConn>::GrainRndPackerPtrVector &
    GougeConfig<TGPckr,TPPckr,TConn>::getGougeGeneratorVector() const
    {
      return m_gougeGenPtrVector;
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    typename GougeConfig<TGPckr,TPPckr,TConn>::GrainRndPackerPtrVector &
    GougeConfig<TGPckr,TPPckr,TConn>::getGougeGeneratorVector()
    {
      return m_gougeGenPtrVector;
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    const typename GougeConfig<TGPckr,TPPckr,TConn>::GeneratorPtrVector &
    GougeConfig<TGPckr,TPPckr,TConn>::getFaultGeneratorVector() const
    {
      return m_faultGenPtrVector;
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    bool GougeConfig<TGPckr,TPPckr,TConn>::areInDifferentFaultBlocks(
      const Particle &p1,
      const Particle &p2
    ) const
    {
      const GeneratorPtrVector &generators = getFaultGeneratorVector();
      if (generators.size() == 2) {
        return
          (
            (generators[0]->contains(p1) && generators[1]->contains(p2))
            ||
            (generators[0]->contains(p2) && generators[1]->contains(p1))
          );
      }
      else if (generators.size() > 2) {
        throw
          std::runtime_error(
            "GougeConfig<TGPckr,TPPckr,TConn>::areInDifferentFaultBlocks: "
            "More than two fault blocks."
          );
      }
      return false;
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    bool GougeConfig<TGPckr,TPPckr,TConn>::isGougeParticle(const Particle &particle) const
    {
      const GrainRndPackerPtrVector &generators = getGougeGeneratorVector();
      for (
        typename GrainRndPackerPtrVector::const_iterator it = generators.begin();
        it != generators.end();
        it++
      )
      {
        if ((*it)->contains(particle)) {
          return true;
        }
      }
      return false;
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    void GougeConfig<TGPckr,TPPckr,TConn>::createConnectionSet()
    {
      //
      // First created connections in the elastic blocks.
      //
      typename NTable::ParticleIterator particleIt = m_nTablePtr->getParticleIterator();
      ConnectionValidator validator = ConnectionValidator(*this, m_prms.getConnectionTolerance());
      while (particleIt.hasNext()) {
        const typename NTable::Particle *pParticle = particleIt.next();
        const typename NTable::ParticleVector neighbours =
          m_nTablePtr->getNeighbourVector(
            pParticle->getPos(),
            pParticle->getRad() + m_prms.getConnectionTolerance()
          );
        for (
          typename NTable::ParticleVector::const_iterator it = neighbours.begin();
          it != neighbours.end();
          it++
        )
        {
          if (validator.isValid(*pParticle, *(*it))) {
            m_connectionSet.insert(
              typename ConnectionSet::value_type(
                pParticle->getID(),
                (*it)->getID(),
                m_prms.getBlockConnectionTag()
              )
            );
          }
        }
      }
      const int numBlockConns = m_connectionSet.size();
      console.Info()
        << "Created " <<  numBlockConns << " connections in "
        << "bonded blocks.\n";

      //
      // Create connections with grains.
      //
      console.Debug()
        << "Prms BBox: " <<  StringUtil::toString(m_prms.getBBox()) << "\n";
      console.Debug()
        << "NTbl BBox: " <<  StringUtil::toString(m_nTablePtr->getBBox()) << "\n";
      for (
        typename GrainRndPackerPtrVector::iterator packerIt = getGougeGeneratorVector().begin();
        packerIt != getGougeGeneratorVector().end();
        packerIt++
      )
      {
        typename GrainRandomPacker::GrainIterator grainIt = (*packerIt)->getGrainIterator();
        while (grainIt.hasNext())
        {
          ConnectionFinder connFinder =
            ConnectionFinder(
              m_prms.getConnectionTolerance(),
              m_prms.getGougeConnectionTag(),
              m_nTablePtr->getBBox(),
              m_nTablePtr->getPeriodicDimensions()
            );
          Grain &g = grainIt.next();
          connFinder.create(g.getParticleIterator());
          typename ConnectionFinder::Iterator connIt = connFinder.getIterator();
          while (connIt.hasNext())
          {
            m_connectionSet.insert(connIt.next());
          }
          if (connFinder.getNumConnections() == 0)
          {
            console.Info()
              << "Found no connections in grain " << g.getId()
              << ":\n";
            typename Grain::ParticleIterator partIt = g.getParticleIterator();
            while (partIt.hasNext())
            {
              console.Info() << StringUtil::toString(partIt.next()) << "\n";
            }
          }
        }
      }
      console.Info()
        << "Created " <<  m_connectionSet.size()-numBlockConns << " connections in "
        << "gouge region.\n";
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    typename GougeConfig<TGPckr,TPPckr,TConn>::ParticleCollection
    GougeConfig<TGPckr,TPPckr,TConn>::getParticleCollection()
    {
      ParticleCollection pCollection(m_particlePoolPtr);
      for (
        typename GeneratorPtrVector::iterator it = m_genPtrVector.begin();
        it != m_genPtrVector.end();
        it++
      )
      {
        ParticleIterator particleIt = (*it)->getParticleIterator();
        while (particleIt.hasNext()) {
          pCollection.insertRef(particleIt.next());
        }
      }

      return pCollection;
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    typename GougeConfig<TGPckr,TPPckr,TConn>::GrainCollection
    GougeConfig<TGPckr,TPPckr,TConn>::getGrainCollection()
    {
      GrainCollection gCollection(m_particlePoolPtr, m_grainPoolPtr);
      for (
        typename GrainRndPackerPtrVector::iterator packerIt =
          getGougeGeneratorVector().begin();
        packerIt != getGougeGeneratorVector().end();
        packerIt++
      )
      {
        GrainIterator grainIt = (*packerIt)->getGrainIterator();
        while (grainIt.hasNext()) {
          gCollection.insertRef(grainIt.next());
        }
      }

      return gCollection;
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    const typename GougeConfig<TGPckr,TPPckr,TConn>::ConnectionSet &
    GougeConfig<TGPckr,TPPckr,TConn>::getConnectionSet() const
    {
      return m_connectionSet;
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    void GougeConfig<TGPckr,TPPckr,TConn>::write(std::ostream &oStream) const
    {
      Vec3 minPt = m_nTablePtr->getBBox().getMinPt();
      Vec3 maxPt = m_nTablePtr->getBBox().getMaxPt();
      if (fabs(maxPt.Z() - minPt.Z()) < (2*m_prms.getMaxRadius())) {
        minPt.Z() = minPt.Z() - m_prms.getMaxRadius() - m_prms.getTolerance();
        maxPt.Z() = maxPt.Z() + m_prms.getMaxRadius() + m_prms.getTolerance();
      }
      
      const BoundingBox geoBBox(minPt, maxPt + m_prms.getTolerance());

      GeometryInfo info = 
        GeometryInfo(
          1.2,
          geoBBox.getMinPt(),
          geoBBox.getMaxPt(),
          m_prms.getPeriodicDimensions(),
          (m_prms.getBBox().getSizes().Z() <= 0.0)
        );
      info.write(oStream);

      /*
       * Some ugliness to ensure that duplicated particles (because of
       * circular boundary conditions) get the correct tag. First, create
       * a set of particles which have already been tagged.
       */
      typename NTable::ParticleIterator particleIt = m_nTablePtr->getParticleIterator();
      typedef std::set<Particle *, IdCompare> ParticleSet;
      ParticleSet taggedParticleSet;
      while (particleIt.hasNext()) {
        taggedParticleSet.insert(particleIt.next());
      }

      /*
       * Now eliminate any particles which were generated with their centre-point
       * lying outside the geo bounding box. Also set the tag in case the particle
       * was a duplicate, created due to circular boundary.
       */
      particleIt = m_nTablePtr->getParticleIterator();
      ParticleSet particleSet;
      while (particleIt.hasNext()) {
        Particle *pParticle = particleIt.next();
        if (geoBBox.contains(pParticle->getPos())) {
          pParticle->setTag((*(taggedParticleSet.find(pParticle)))->getTag());
          particleSet.insert(pParticle);
        }
      }
      
      /*
       * Write particles to the stream.
       */
      oStream 
        << "\n" 
        << "BeginParticles"
        << "\n"
        << "Simple"
        << "\n"
        << particleSet.size()
        << "\n";
      const int precision = 12;
      GeoParticleWriter particleVisitor(oStream, precision);
      for (
        typename ParticleSet::const_iterator it = particleSet.begin();
        it != particleSet.end();
        it++
      )
      {
        particleVisitor.visitParticle(*(*it));
      }

      oStream << "EndParticles\n" << "BeginConnect\n";
      oStream << getConnectionSet().size() << "\n";
      oStream.flush();

      GeoConnectionWriter connectionVisitor(oStream);
      visitConnections(connectionVisitor);
      oStream << "EndConnect";
      oStream.flush();
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    void GougeConfig<TGPckr,TPPckr,TConn>::tagGougeParticles(int tag)
    {
      for (
        typename GrainRndPackerPtrVector::iterator it = m_gougeGenPtrVector.begin();
        it != m_gougeGenPtrVector.end();
        it++
      )
      {
        typename GrainRandomPacker::ParticleIterator particleIt =
          (*it)->getParticleIterator();
        while (particleIt.hasNext()) {
          particleIt.next().setTag(tag);
        }
      }
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    void GougeConfig<TGPckr,TPPckr,TConn>::tagRndBlockParticles(int tag)
    {
      for (
        typename GeneratorPtrVector::iterator it = m_faultGenPtrVector.begin();
        it != m_faultGenPtrVector.end();
        it++
      )
      {
        ParticleIterator particleIt = (*it)->getParticleIterator();
        while (particleIt.hasNext()) {
          particleIt.next().setTag(tag);
        }
      }
    }

    template <typename TGPckr,typename TPPckr,typename TConn>
    void GougeConfig<TGPckr,TPPckr,TConn>::tagDrivingPlateParticles(
      int lowDrivingTag,
      int highDrivingTag,
      double distanceFromBBoxEdge
    )
    {
      ParticleCollection particleCollection = getParticleCollection();
      const BoundingBox bBox = particleCollection.getParticleBBox();
      
      const int    idx     = this->m_prms.getOrientationIndex();
      const double maxLow  = bBox.getMinPt()[idx] + distanceFromBBoxEdge;
      const double minHigh = bBox.getMaxPt()[idx] - distanceFromBBoxEdge;

      int lowTagCount = 0;
      int highTagCount = 0;
      typename ParticleCollection::ParticleIterator particleIt =
        particleCollection.getParticleIterator();
      while (particleIt.hasNext())
      {
        Particle &particle = particleIt.next();
        const double dimPos = particle.getPos()[idx];
        const double radius = particle.getRad();
        
        if (dimPos - radius <= maxLow) {
          particle.setTag(lowDrivingTag);
          lowTagCount++;
        }
        if (dimPos + radius >= minHigh) {
          particle.setTag(highDrivingTag);
          highTagCount++;
        }
      }
      console.Info() << "Tagged " << lowTagCount << " particles with " << lowDrivingTag << "\n";
      console.Info() << "Tagged " << highTagCount << " particles with " << highDrivingTag << "\n";
    }
  }
}
