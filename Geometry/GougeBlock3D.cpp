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
#include "Geometry/GougeBlock3D.h"
#include "Geometry/GeometryInfo.h"
#include "Geometry/CircularNeighbourTable.h"
#include "Geometry/RegularBlockGenerator.h"
#include "Geometry/RandomBlockGenerator.h"
#include "Geometry/GridIterator.h"
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

    double PackingInfo::getMinRadius() const
    {
      return m_minRadius;
    }

    double PackingInfo::getMaxRadius() const
    {
      return m_maxRadius;
    }

    const BoolVector &PackingInfo::getPeriodicDimensions() const
    {
      return m_periodicDimensions;
    }

    //==========================================================================

    ParticleBlockPrms::ParticleBlockPrms()
      : m_size(0.0),
        m_minRadius(0.0),
        m_maxRadius(0.0)
    {
    }

    ParticleBlockPrms::ParticleBlockPrms(
      double size,
      double minRadius,
      double maxRadius
    ) : m_size(size),
        m_minRadius(minRadius),
        m_maxRadius(maxRadius)
    {
    }

    ParticleBlockPrms::~ParticleBlockPrms()
    {
    }
    
    //==========================================================================

    GougeBlockPrms::GougeBlockPrms()
      : m_bBox(Vec3::ZERO, Vec3::ZERO),
        m_padRadius(0.0),
        m_orientation(XZ),
        m_faultPrms(),
        m_gougePrms(),
        m_periodicDimensions(3, false),
        m_maxInsertionFailures(50),
        m_tolerance(DBL_EPSILON*128),
        m_connectionTolerance(DBL_EPSILON*128*10)
    {
    }

    GougeBlockPrms::GougeBlockPrms(
      const BoundingBox       &bBox,
      double                  padRadius,
      Orientation             orientation,
      const ParticleBlockPrms &faultRegionPrms,
      const ParticleBlockPrms &gougeRegionPrms,
      const BoolVector        &periodicDimensions,
      int                     maxInsertionFailures,
      double                  tolerance ,
      double                  connectionTolerance
    ) : m_bBox(bBox),
        m_padRadius(padRadius),
        m_orientation(orientation),
        m_faultPrms(faultRegionPrms),
        m_gougePrms(gougeRegionPrms),
        m_periodicDimensions(periodicDimensions),
        m_maxInsertionFailures(maxInsertionFailures),
        m_tolerance(tolerance),
        m_connectionTolerance(connectionTolerance)
    {
      m_bBox = GridIterator(bBox, getMaxRadius()).getSphereBBox();
    }

    GougeBlockPrms::~GougeBlockPrms()
    {
    }

    const BoundingBox &GougeBlockPrms::getBBox() const
    {
      return m_bBox;
    }

    int GougeBlockPrms::getMaxInsertionFailures() const
    {
      return m_maxInsertionFailures;
    }
    
    const BoolVector &GougeBlockPrms::getPeriodicDimensions() const
    {
      return m_periodicDimensions;
    }

    Orientation GougeBlockPrms::getOrientation() const
    {
      return m_orientation;
    }
    
    double GougeBlockPrms::getTolerance() const
    {
      return m_tolerance;
    }

    double GougeBlockPrms::getConnectionTolerance() const
    {
      return m_connectionTolerance;
    }
    
    int GougeBlockPrms::getOrientationIndex() const
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

    /*!
      generate a new bounding box on one side of the centre between two given distances

      \param d1 near distance
      \param d2 far distance
    */
    BoundingBox GougeBlockPrms::cutFromCentre(double d1, double d2) const
    {
      const int idx = getOrientationIndex();
      const BoundingBox bBox = getBBox(); 

      // get near/far points 
      const double centerPoint=0.5*(bBox.getMaxPt()[idx] + bBox.getMinPt()[idx]);
      const double cmp1 = d1 + centerPoint;
      const double cmp2 = d2 + centerPoint;
      
      // min/max point
      Vec3 minPt = bBox.getMinPt();
      Vec3 maxPt = bBox.getMaxPt();

      minPt[idx] = std::min(cmp1, cmp2);
      maxPt[idx] = std::max(cmp1, cmp2);

      // get maximum regular grid for box of size maxPt-minPt
      const BoundingBox gridBBox = GridIterator(BoundingBox(Vec3(0.0,0.0,0.0), maxPt-minPt), getMaxRadius()).getSphereBBox();
      const double oldDim= maxPt[idx]-minPt[idx];
      const double newDim= gridBBox.getSizes()[idx];
      const double diffDim=oldDim-newDim;
      // adjust near side
      if(d1>0.0){ // "positive" side -> adjust minPt
	minPt[idx]+=diffDim;
      } else { // "negative" side -> adjust maxPt
	maxPt[idx]-=diffDim;
      }
      BoundingBox returnBBox = BoundingBox(minPt, maxPt);

      return returnBBox;
    }

    double GougeBlockPrms::getRegularBlockRadius() const
    {
      return m_padRadius;
    }

    double GougeBlockPrms::getFaultMinRadius() const
    {
      return m_faultPrms.m_minRadius;
    }

    double GougeBlockPrms::getFaultMaxRadius() const
    {
      return m_faultPrms.m_maxRadius;
    }

    double GougeBlockPrms::getGougeMinRadius() const
    {
      return m_gougePrms.m_minRadius;
    }

    double GougeBlockPrms::getGougeMaxRadius() const
    {
      return m_gougePrms.m_maxRadius;
    }

    double GougeBlockPrms::getOrientationSize() const
    {
      return 
        (
          getBBox().getMaxPt()[getOrientationIndex()]
          -
          getBBox().getMinPt()[getOrientationIndex()]
        );
    }

    BoundingBoxVector GougeBlockPrms::getRegularBBoxVector() const
    {
      BoundingBoxVector bBoxVector;
      if (m_gougePrms.m_size + 2*m_faultPrms.m_size < getOrientationSize()) {
        bBoxVector.reserve(2);
        bBoxVector.push_back(
          cutFromCentre(
            -(m_gougePrms.m_size/2.0 + m_faultPrms.m_size),
            -getOrientationSize()/2.0
          )
        );
        bBoxVector.push_back(
          cutFromCentre(
            m_gougePrms.m_size/2.0 + m_faultPrms.m_size,
            getOrientationSize()/2.0
          )
        );
      }
      return bBoxVector;
    }

    PackingInfoVector GougeBlockPrms::getGougePackingInfoVector() const
    {
      PackingInfoVector packingInfoVector;
      if (m_gougePrms.m_size > 0.0) {
        Vec3 overlap = Vec3::ZERO;
        overlap[getOrientationIndex()] = m_faultPrms.m_maxRadius;
        BoundingBox bBox =
          cutFromCentre(
             m_gougePrms.m_size/2.0,
            -m_gougePrms.m_size/2.0
          );
        
        packingInfoVector.push_back(
          PackingInfo(
            BoundingBox(bBox.getMinPt() - overlap, bBox.getMaxPt() + overlap),
            getPeriodicDimensions(),
            getOrientation(),
            getGougeMinRadius(),
            getGougeMaxRadius()
          )
        );
      }
      return packingInfoVector;
    }

    PackingInfoVector GougeBlockPrms::getFaultPackingInfoVector() const
    {
      PackingInfoVector packingInfoVector;
      if (m_faultPrms.m_size > 0.0) {
        packingInfoVector.reserve(2);
        Vec3 overlap = Vec3::ZERO;
        overlap[getOrientationIndex()] = m_padRadius;
        BoundingBox bBox = 
          cutFromCentre(
            -m_gougePrms.m_size/2.0,
            -(m_gougePrms.m_size/2.0 + m_faultPrms.m_size)
          );
        packingInfoVector.push_back(
          PackingInfo(
            BoundingBox(bBox.getMinPt() - overlap, bBox.getMaxPt()),
            getPeriodicDimensions(),
            getOrientation(),
            getFaultMinRadius(),
            getFaultMaxRadius()
          )
        );
        bBox =
          cutFromCentre(
            m_gougePrms.m_size/2.0,
            m_gougePrms.m_size/2.0 + m_faultPrms.m_size
          );

        packingInfoVector.push_back(
          PackingInfo(
            BoundingBox(bBox.getMinPt(), bBox.getMaxPt() + overlap),
            getPeriodicDimensions(),
            getOrientation(),
            getFaultMinRadius(),
            getFaultMaxRadius()
          )
        );
      }
      return packingInfoVector;
    }

    double GougeBlockPrms::getMaxRadius() const
    {
      return 
        std::max(
          m_padRadius,
          std::max(m_faultPrms.m_maxRadius, m_gougePrms.m_maxRadius)
        );
    }

    double GougeBlockPrms::getMinRadius() const
    {
      return 
        std::min(
          m_padRadius,
          std::min(m_faultPrms.m_minRadius, m_gougePrms.m_minRadius)
        );
    }

    bool GougeBlockPrms::is2d() const
    {
      return (getBBox().getSizes()[2] == 0.0);
    }

    //==========================================================================
    
    int GougeBlock3D::getNumParticles() const
    {
      int numParticles = 0;
      for (
        GeneratorPtrVector::const_iterator it = m_genPtrVector.begin();
        it != m_genPtrVector.end();
        it++
      )
      {
        numParticles += (*it)->getNumParticles();
      }
      return numParticles;
    }

    GougeBlock3D::GougeBlock3D(const GougeBlockPrms &prms)
      : m_nTablePtr(),
        m_prms(prms),
        m_interactionSet(),
        m_gougeGenPtrVector(),
        m_genPtrVector(),
        m_particlePoolPtr(new ParticlePool(2048)),
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

    void GougeBlock3D::createRegularBlockGenerators()
    {      
      BoundingBoxVector bBoxVector = m_prms.getRegularBBoxVector();
      for (
        BoundingBoxVector::const_iterator it = bBoxVector.begin();
        it != bBoxVector.end();
        it++
      ) {
	std::cout << "regular block bounding box is :" << *it << std::endl; 
        GeneratorPtr genPtr = 
          GeneratorPtr(
            new RegularBlockGenerator(
              *(m_nTablePtr.get()),
              *(m_particlePoolPtr.get()),
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

    void GougeBlock3D::createFaultBlockGenerators()
    {      
      PackingInfoVector packingInfoVector = m_prms.getFaultPackingInfoVector();
      for (
        PackingInfoVector::const_iterator it = packingInfoVector.begin();
        it != packingInfoVector.end();
        it++
      ) {
        GeneratorPtr genPtr = 
          GeneratorPtr(
            new RandomBlockGenerator(
              *(m_nTablePtr.get()),
              *(m_particlePoolPtr.get()),
              it->getBBox(),
              it->getPeriodicDimensions(),
              m_prms.getTolerance(),
              it->getMinRadius(),
              it->getMaxRadius(),
              it->getFitPlaneVector(),
              m_prms.getMaxInsertionFailures()
            )
          );
        m_genPtrVector.push_back(genPtr);
        m_faultGenPtrVector.push_back(genPtr);
      }
    }

    void GougeBlock3D::createGougeBlockGenerators()
    {
      PackingInfoVector packingInfoVector = m_prms.getGougePackingInfoVector();
      for (
        PackingInfoVector::const_iterator it = packingInfoVector.begin();
        it != packingInfoVector.end();
        it++
      ) {
        GeneratorPtr genPtr = 
          GeneratorPtr(
            new RandomBlockGenerator(
              *(m_nTablePtr.get()),
              *(m_particlePoolPtr.get()),
              it->getBBox(),
              it->getPeriodicDimensions(),
              m_prms.getTolerance(),
              it->getMinRadius(),
              it->getMaxRadius(),
              it->getFitPlaneVector(),
              m_prms.getMaxInsertionFailures()
            )
          );
        m_genPtrVector.push_back(genPtr);
        m_gougeGenPtrVector.push_back(genPtr);
      }
    }

    GougeBlock3D::~GougeBlock3D()
    {
    }

    void GougeBlock3D::generate()
    {
      // setup block generators
      createRegularBlockGenerators();
      createFaultBlockGenerators();
      createGougeBlockGenerators();
 
      // use block generators  
      console.Info() << "bbox = " << m_prms.getBBox().getMinPt() << " " << m_prms.getBBox().getMaxPt() << "\n";
      for (
        GeneratorPtrVector::iterator it = m_genPtrVector.begin();
        it != m_genPtrVector.end();
        it++
      )
      {
        (*it)->generate();
      }
      createInteractionSet();
    }

    void GougeBlock3D::writeToFile(const std::string &fileName) const
    {
      std::ofstream fStream(fileName.c_str());
      write(fStream);
    }

    const GougeBlock3D::GeneratorPtrVector &GougeBlock3D::getGougeGeneratorVector() const
    {
      return m_gougeGenPtrVector;
    }

    const GougeBlock3D::GeneratorPtrVector &GougeBlock3D::getFaultGeneratorVector() const
    {
      return m_faultGenPtrVector;
    }

    bool GougeBlock3D::areInDifferentFaultBlocks(
      const SimpleParticle &p1,
      const SimpleParticle &p2
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
        throw std::runtime_error("GougeBlock3D::areInDifferentFaultBlocks: More than two fault blocks.");
      }
      return false;
    }

    bool GougeBlock3D::isGougeParticle(const SimpleParticle &particle) const
    {
      const GeneratorPtrVector &generators = getGougeGeneratorVector();
      for (
        GeneratorPtrVector::const_iterator it = generators.begin();
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

    class InteractionValidator
    {
    public:
      inline InteractionValidator(const GougeBlock3D &gougeBlock, double tolerance)
        : m_pGougeBlock(&gougeBlock),
          m_tolerance(tolerance)
      {
      }

      inline bool isValid(const SimpleParticle &p1, const SimpleParticle &p2) const
      {
        return 
          (
            (p1.getID() < p2.getID())
            &&
            ((p1.getPos() - p2.getPos()).norm() < (m_tolerance + (p1.getRad() + p2.getRad())))
            &&
            ((!m_pGougeBlock->isGougeParticle(p1)) && (!m_pGougeBlock->isGougeParticle(p2)))
            &&
            ((!m_pGougeBlock->areInDifferentFaultBlocks(p1, p2)))
          );
      }

    private:
      const GougeBlock3D *m_pGougeBlock;
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
      
      void visitSimpleParticle(const SimpleParticle &particle) const
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

    class GeoInteractionWriter
    {
    public:
      GeoInteractionWriter(std::ostream &oStream)
        : m_pOStream(&oStream)
      {
      }

      void visitBasicInteraction(const BasicInteraction &interaction)
      {
        (*m_pOStream)
          << interaction.first()  << " "
          << interaction.second() << " "
          << 0 << "\n";
      }

    private:
      std::ostream *m_pOStream;
      int          m_precision;
    };

    void GougeBlock3D::createInteractionSet()
    {
      NTable::ParticleIterator particleIt = m_nTablePtr->getParticleIterator();
      InteractionValidator validator = InteractionValidator(*this, m_prms.getConnectionTolerance());
      while (particleIt.hasNext()) {
        const NTable::Particle *pParticle = particleIt.next();
        const NTable::ParticleVector neighbours =
          m_nTablePtr->getNeighbourVector(
            pParticle->getPos(),
            pParticle->getRad() + m_prms.getConnectionTolerance()
          );
        for (
          NTable::ParticleVector::const_iterator it = neighbours.begin();
          it != neighbours.end();
          it++
        )
        {
          if (validator.isValid(*pParticle, *(*it))) {
            m_interactionSet.insert(
              BasicInteraction(pParticle->getID(), (*it)->getID())
            );
          }
        }
      }
    }

    const GougeBlock3D::InteractionSet &GougeBlock3D::getInteractionSet() const
    {
      return m_interactionSet;
    }

    class IdCompare
    {
    public:
      bool operator()(const SimpleParticle *p1, const SimpleParticle *p2) const
      {
        return (p1->getID() < p2->getID());
      }
    };
    
    void GougeBlock3D::write(std::ostream &oStream) const
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
      NTable::ParticleIterator particleIt = m_nTablePtr->getParticleIterator();
      typedef std::set<SimpleParticle *, IdCompare> ParticleSet;
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
        SimpleParticle *pParticle = particleIt.next();
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
      for (ParticleSet::const_iterator it = particleSet.begin(); it != particleSet.end(); it++)
      {
        particleVisitor.visitSimpleParticle(*(*it));
      }

      oStream << "EndParticles\n" << "BeginConnect\n";
      oStream << getInteractionSet().size() << "\n";      
      oStream.flush();

      GeoInteractionWriter interactionVisitor(oStream);
      visitInteractions(interactionVisitor);
      oStream << "EndConnect";
      oStream.flush();
    }

    void GougeBlock3D::tagGougeParticles(int tag)
    {
      for (
        GeneratorPtrVector::iterator it = m_gougeGenPtrVector.begin();
        it != m_gougeGenPtrVector.end();
        it++
      )
      {
        BlockGenerator::ParticleIterator particleIt = (*it)->getParticleIterator();
        while (particleIt.hasNext()) {
          particleIt.next()->setTag(tag);
        }
      }
    }

    void GougeBlock3D::tagFaultParticles(int tag)
    {
      for (
        GeneratorPtrVector::iterator it = m_faultGenPtrVector.begin();
        it != m_faultGenPtrVector.end();
        it++
      )
      {
        BlockGenerator::ParticleIterator particleIt = (*it)->getParticleIterator();
        while (particleIt.hasNext()) {
          particleIt.next()->setTag(tag);
        }
      }
    }
    
    void GougeBlock3D::tagDrivingPlateParticles(
      int lowDrivingTag,
      int highDrivingTag,
      double distanceFromBBoxEdge
    )
    {
      const int    idx     = this->m_prms.getOrientationIndex();
      const double maxLow  = m_prms.getBBox().getMinPt()[idx] + distanceFromBBoxEdge;
      const double minHigh = m_prms.getBBox().getMaxPt()[idx] - distanceFromBBoxEdge;

      int lowTagCount = 0;
      int highTagCount = 0;
      for (
        GeneratorPtrVector::iterator it = m_regularGenPtrVector.begin();
        it != m_regularGenPtrVector.end();
        it++
      )
      {
        console.Info() << (*it)->getBBox().getMinPt() << " " << (*it)->getBBox().getMaxPt() << "\n";
        BlockGenerator::ParticleIterator particleIt = (*it)->getParticleIterator();
        while (particleIt.hasNext()) {
          SimpleParticle *pParticle = particleIt.next();
          const double dimPos = pParticle->getPos()[idx];
          const double radius = pParticle->getRad();
          
          if (dimPos - radius <= maxLow) {
            pParticle->setTag(lowDrivingTag);
            lowTagCount++;
          }
          if (dimPos + radius >= minHigh) {
            pParticle->setTag(highDrivingTag);
            highTagCount++;
          }
        }
      }
      console.Info() << "Tagged " << lowTagCount << " particles with " << lowDrivingTag << "\n";
      console.Info() << "Tagged " << highTagCount << " particles with " << highDrivingTag << "\n";
    }
  }
}
