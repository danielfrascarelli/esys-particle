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


#include "Geometry/GeometryInfo.h"

// --- STL includes --- 
#include <stdexcept>
#include <sstream>
#include <algorithm> 

namespace esys
{
  namespace lsm
  {
    class GeometryInfo::Impl
    {
    public:
      Impl();
      
      Impl(
        float            version,
        const Vec3       &bBoxMin,
        const Vec3       &bBoxMax,
        const BoolVector &periodicDimensions,
        bool             is2d = false
      );

      Impl &operator=(const Impl &impl);

      bool operator==(const Impl &impl) const;

      ~Impl();

      void read(std::istream &iStream);

      void write(std::ostream &oStream) const;
      void writeWithoutVersion(std::ostream &oStream) const;

      float      m_version;
      Vec3       m_bBoxMin;
      Vec3       m_bBoxMax;
      /**
       * Indicates the dimensions which are periodic, element 0 indicates
       * whether the x-dim is periodic, element 1 indicates whether the y-dim
       * is periodic and element 2 indicates whether the z-dim is periodic.
       */
      BoolVector m_periodicDimensions;
      bool       m_is2d;      
    };

    GeometryInfo::Impl::Impl()
      : m_version(0.0),
        m_bBoxMin(),
        m_bBoxMax(),
        m_periodicDimensions(3, false),
        m_is2d(false)        
    {
    }

    GeometryInfo::Impl::Impl(
      float            version,
      const Vec3       &bBoxMin,
      const Vec3       &bBoxMax,
      const BoolVector &periodicDimensions,
      bool             is2d
    )
      : m_version(version),
        m_bBoxMin(bBoxMin),
        m_bBoxMax(bBoxMax),
        m_periodicDimensions(periodicDimensions),
        m_is2d(is2d)
    {
    }

    GeometryInfo::Impl::~Impl()
    {
    }

    GeometryInfo::Impl &GeometryInfo::Impl::operator=(const GeometryInfo::Impl &impl)
    {
      m_version            = impl.m_version;
      m_bBoxMin            = impl.m_bBoxMin;
      m_bBoxMax            = impl.m_bBoxMax;
      m_periodicDimensions = impl.m_periodicDimensions;
      m_is2d               = impl.m_is2d;      

      return *this;
    }

    bool GeometryInfo::Impl::operator==(const GeometryInfo::Impl &impl) const
    {
      return
        (
          (m_version == impl.m_version)
          &&
          (m_bBoxMin == impl.m_bBoxMin)
          &&
          (m_bBoxMax == impl.m_bBoxMax)
          &&
          (m_is2d == impl.m_is2d)
          &&
          (m_periodicDimensions == impl.m_periodicDimensions)
        );
    }

    /*!
      Write complete geometry info to output stream

      \param oStream the output stream
    */
    void GeometryInfo::Impl::write(std::ostream &oStream) const
    {
      oStream << "LSMGeometry " << m_version << "\n";
      writeWithoutVersion(oStream);
    }

    /*!
      Write geometry info to output stream, except the version nr. 
      of the .geo file. This output can't be read by the "read" function.
      It is used for the snapshot headers.

      \param oStream the output stream
    */
    void GeometryInfo::Impl::writeWithoutVersion(std::ostream &oStream) const
    {
      oStream << "BoundingBox ";
      oStream 
        << m_bBoxMin.X() << " "
        << m_bBoxMin.Y() << " "
        << m_bBoxMin.Z() << " ";
      oStream 
        << m_bBoxMax.X() << " "
        << m_bBoxMax.Y() << " "
        << m_bBoxMax.Z() << "\n";

      oStream << "PeriodicBoundaries ";
      oStream 
        << m_periodicDimensions[0] << " "
        << m_periodicDimensions[1] << " "
        << m_periodicDimensions[2];

      oStream
	<< "\nDimension "
	<< (m_is2d ? "2D" : "3D");
    }

    void GeometryInfo::Impl::read(std::istream &iStream)
    {
      Impl impl;

      // read and check file type marker
      std::string fileType;
      iStream >> fileType;
      size_t ftl=std::min(fileType.size(),size_t(4));
      std::string fileTypeHead=fileType.substr(0,ftl);
      if ((fileType != "LSMGeometry") && (fileTypeHead!="ESyS")) { // wrong file type -> bail out
        throw runtime_error("Unrecognised file type " + fileType + " expected LSMGeometry or ESyS-Particle.");
      }

      // read and check version number
      iStream >> impl.m_version;
      if((fileType=="LSMGeometry") && ((impl.m_version != 1.1f) && (impl.m_version != 1.2f) && (impl.m_version != 1.3f))){
        std::stringstream msg;
        msg
          << "Can only read version 1.1, 1.2 or 1.3 geometry files, this is version "
          << impl.m_version;
        throw std::runtime_error(msg.str().c_str());
      } else if (fileTypeHead=="ESyS") {
	impl.m_version=1.2f; // temporary hack 
      }

      // read boundary box
      string bbx;
      iStream >> bbx;
      // check boundary box marker
      if (bbx != "BoundingBox") {
        throw std::runtime_error("Expected BoundingBox, got " + bbx);;
      }
      iStream >> impl.m_bBoxMin.X() >> impl.m_bBoxMin.Y() >> impl.m_bBoxMin.Z();
      iStream >> impl.m_bBoxMax.X() >> impl.m_bBoxMax.Y() >> impl.m_bBoxMax.Z();

      // read boundary periodicity
      string bdry;
      iStream >> bdry;
      // check boundary periodicity marker
      if(bdry != "PeriodicBoundaries")
      {
        throw std::runtime_error("Expected PeriodicBoundaries, got " + bdry);
      }
      for (int i = 0; i < 3; i++) {
        bool isPeriodic = false;
        iStream >> isPeriodic;
        impl.m_periodicDimensions[i] = isPeriodic;
      }

      // if 1.2, read in 2D or 3D
      if(impl.m_version >= 1.2f){
        std::string dims;
        iStream >> dims;
        if(dims != "Dimension") {
          throw std::runtime_error("Expected 'Dimension', got '" + dims + "'");
        }

        std::string r2d;        
        iStream >> r2d;
        if((r2d == "2D") || (r2d == "2d")){
          impl.m_is2d = true;
        } else {
          impl.m_is2d = false;
        }
      } else {
        impl.m_is2d = true;
      }

	  // if 1.3, read generation version info -> currently not used
	  if(impl.m_version == 1.3f){
        std::string ggv;
        iStream >> ggv;
        if(ggv != "GengeoVersion") {
          throw std::runtime_error("Expected 'GengeoVersion', got '" + ggv + "'");
        }

        std::string version, bzr_rev;        
        iStream >> version >> bzr_rev;
      } 
      
      *this = impl;
    }

    GeometryInfo::GeometryInfo()
      : m_pImpl(new GeometryInfo::Impl())
    {
    }

    GeometryInfo::GeometryInfo(
      float            version,
      const Vec3       &bBoxMin,
      const Vec3       &bBoxMax,
      const BoolVector &periodicDimensions,
      bool             is2d
    )
      : m_pImpl(new Impl(version, bBoxMin, bBoxMax, periodicDimensions, is2d))
    {
    }

    GeometryInfo::GeometryInfo(const GeometryInfo &geoInfo)
      : m_pImpl(new GeometryInfo::Impl(*(geoInfo.m_pImpl)))
    {
    }

    GeometryInfo &GeometryInfo::operator=(const GeometryInfo &geoInfo)
    {
      *m_pImpl = *geoInfo.m_pImpl;
      return *this;
    }

    GeometryInfo::~GeometryInfo()
    {
      delete m_pImpl;
    }

    bool GeometryInfo::operator==(const GeometryInfo &geoInfo) const
    {
      return (*m_pImpl == *geoInfo.m_pImpl);
    }

    bool GeometryInfo::hasAnyPeriodicDimensions() const
    {
      return 
        (
          ((m_pImpl->m_periodicDimensions.size() > 0) && (m_pImpl->m_periodicDimensions[0]))
          ||
          ((m_pImpl->m_periodicDimensions.size() > 1) && (m_pImpl->m_periodicDimensions[1]))
          ||
          ((m_pImpl->m_periodicDimensions.size() > 2) && (m_pImpl->m_periodicDimensions[2]))
        );
    }

    bool GeometryInfo::is2d() const
    {
      return m_pImpl->m_is2d;
    }

    void GeometryInfo::set_is2d(bool do2d)
    {
      m_pImpl->m_is2d = do2d;
    }

    Vec3Vector GeometryInfo::getBBoxCorners() const
    {
      Vec3Vector corners;
      corners.push_back(m_pImpl->m_bBoxMin);
      corners.push_back(m_pImpl->m_bBoxMax);

      return corners;
    }

    Vec3 GeometryInfo::getMinBBoxCorner() const
    {
      return m_pImpl->m_bBoxMin;
    }

    Vec3 GeometryInfo::getMaxBBoxCorner() const
    {
      return m_pImpl->m_bBoxMax;
    }

    IntVector GeometryInfo::getPeriodicDimensions() const
    {
      return 
        IntVector(
          m_pImpl->m_periodicDimensions.begin(),
          m_pImpl->m_periodicDimensions.end()
        );
    }

    void GeometryInfo::setPeriodicDimensions(BoolVector periodicDimensions)
    {
      m_pImpl->m_periodicDimensions = periodicDimensions;
    }

    void GeometryInfo::setLsmGeoVersion(float version)
    {
      m_pImpl->m_version = version;
    }
      
    /*!
	get the version of the loaded geometry
    */    
    float GeometryInfo::getLsmGeoVersion() const
    {
      return m_pImpl->m_version;
    }
    
    void GeometryInfo::read(std::istream &iStream)
    {
      m_pImpl->read(iStream);
    }

    void GeometryInfo::write(std::ostream &oStream) const
    {
      m_pImpl->write(oStream);
    }

    void GeometryInfo::writeWithoutVersion(std::ostream &oStream) const
    {
      m_pImpl->writeWithoutVersion(oStream);
    }

    void GeometryInfo::setBBox(const Vec3 &min, const Vec3 &max)
    {
      m_pImpl->m_bBoxMin = min;
      m_pImpl->m_bBoxMax = max;
    }

	/*!
	 check if a model with the geometry described in the GeoInfo given in the argument
	 can be loaded into a pre-existing model 
	
	 \param gi the new geometry info
	 */
	bool GeometryInfo::isCompatible(const GeometryInfo& gi) const
	{
		bool res=true;
		
		// check if circular dimensions agree
		res=res && (m_pImpl->m_periodicDimensions[0]==gi.m_pImpl->m_periodicDimensions[0]) 
			&& (m_pImpl->m_periodicDimensions[1]==gi.m_pImpl->m_periodicDimensions[1]) 
			&& (m_pImpl->m_periodicDimensions[2]==gi.m_pImpl->m_periodicDimensions[2]); 
		// check if min/max in circular dimensions agree and 
		// if new (argument) bounding box fits into old bbx in the non-circular dimensions
		// x
		if(m_pImpl->m_periodicDimensions[0]) {
			res = res && (m_pImpl->m_bBoxMin[0]==gi.m_pImpl->m_bBoxMin[0])
				&& (m_pImpl->m_bBoxMax[0]==gi.m_pImpl->m_bBoxMax[0]);
		} else {
			res = res && (m_pImpl->m_bBoxMin[0]<=gi.m_pImpl->m_bBoxMin[0])
				&& (m_pImpl->m_bBoxMax[0]>=gi.m_pImpl->m_bBoxMax[0]);
		}
		// y
 		if(m_pImpl->m_periodicDimensions[1]) {
			res = res && (m_pImpl->m_bBoxMin[1]==gi.m_pImpl->m_bBoxMin[1])
				&& (m_pImpl->m_bBoxMax[1]==gi.m_pImpl->m_bBoxMax[1]);
		} else {
			res = res && (m_pImpl->m_bBoxMin[1]<=gi.m_pImpl->m_bBoxMin[1])
				&& (m_pImpl->m_bBoxMax[1]>=gi.m_pImpl->m_bBoxMax[1]);
		}
		// z
		if(m_pImpl->m_periodicDimensions[2]) {
			res = res && (m_pImpl->m_bBoxMin[2]==gi.m_pImpl->m_bBoxMin[2])
				&& (m_pImpl->m_bBoxMax[2]==gi.m_pImpl->m_bBoxMax[2]);
		} else {
			res = res && (m_pImpl->m_bBoxMin[2]<=gi.m_pImpl->m_bBoxMin[2])
				&& (m_pImpl->m_bBoxMax[2]>=gi.m_pImpl->m_bBoxMax[2]);
		}
		
		return res;
	}

    
	/*!
		Check if the geometrical dimensions in two GeometryInfo objects are identical
	
		\param gi the new geometry info
	 */
	bool GeometryInfo::isIdenticalGeometry(const GeometryInfo& gi) const
	{
		bool res=true;
		
		// check if circular dimensions agree
		res=res && (m_pImpl->m_periodicDimensions[0]==gi.m_pImpl->m_periodicDimensions[0]) 
			&& (m_pImpl->m_periodicDimensions[1]==gi.m_pImpl->m_periodicDimensions[1]) 
			&& (m_pImpl->m_periodicDimensions[2]==gi.m_pImpl->m_periodicDimensions[2]); 
		// check if dimensions agree 
		// x
		res = res && (m_pImpl->m_bBoxMin[0]==gi.m_pImpl->m_bBoxMin[0])
				&& (m_pImpl->m_bBoxMax[0]==gi.m_pImpl->m_bBoxMax[0]);
		// y
 		res = res && (m_pImpl->m_bBoxMin[1]==gi.m_pImpl->m_bBoxMin[1])
				&& (m_pImpl->m_bBoxMax[1]==gi.m_pImpl->m_bBoxMax[1]);
		// z
		res = res && (m_pImpl->m_bBoxMin[2]==gi.m_pImpl->m_bBoxMin[2])
				&& (m_pImpl->m_bBoxMax[2]==gi.m_pImpl->m_bBoxMax[2]);
		
		return res;
	}

	
    std::ostream &operator<<(std::ostream &oStream, const GeometryInfo &geoInfo)
    {
      geoInfo.write(oStream);
      return oStream;
    }

    std::istream &operator<<(std::istream &iStream, GeometryInfo &geoInfo)
    {
      geoInfo.read(iStream);
      return iStream;
    }
  }
}
