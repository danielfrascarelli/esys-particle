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


#ifndef ESYS_LSMGEOMETRYINFO_H
#define ESYS_LSMGEOMETRYINFO_H

#include "Foundation/vec3.h"

#include <vector>
#include <iostream>

namespace esys
{
  namespace lsm
  {
    typedef std::vector<bool> BoolVector;
    typedef std::vector<Vec3> Vec3Vector;
    typedef std::vector<int>  IntVector;

    /**
     * Container class for geometry meta-info.
     */
    class GeometryInfo
    {
    public:
      /**
       *
       */
      GeometryInfo();

      /**
       *
       */
      GeometryInfo(
        float            version,
        const Vec3       &bBoxMin,
        const Vec3       &bBoxMax,
        const BoolVector &periodicDimensions,
        bool             is2d = false
      );

      GeometryInfo(const GeometryInfo &geoInfo);

      GeometryInfo &operator=(const GeometryInfo &geoInfo);

      ~GeometryInfo();

      bool operator==(const GeometryInfo &geoInfo) const;

      /**
       * Sets the bounding box for geometry data.
       */
      void setBBox(const Vec3 &min, const Vec3 &max);

      /**
       * Returns true if any of the x, y or z dimensions
       * have been specified as periodic.
       */
      bool hasAnyPeriodicDimensions() const;

      /**
       * Returns true info indicates two-dimensional particle data.
       */
      bool is2d() const;

      /**
       * Set 2-D information to true if the particle data are two-dimensional; 
       * otherwise set to false.
       */
      void set_is2d(bool do2d);
      
      /**
       * Returns two corner points of bounding box.
       */
      Vec3Vector getBBoxCorners() const;
      Vec3 getMinBBoxCorner() const;
      Vec3 getMaxBBoxCorner() const;

      /**
       * Returns the periodic dimensions.
       */
      IntVector getPeriodicDimensions() const;

      /**
       * Set the periodicity of the x, y and z dimensions.
       */
      void setPeriodicDimensions(BoolVector periodicDimensions);

      /**
       * Set the LSMGeometry version for use in geometry files.
       */
      void setLsmGeoVersion(float version);
      float getLsmGeoVersion() const;
      /**
       * Parses specified istream and assigns to this object.
       */
      void read(std::istream &iStream);

      /**
       * Writes to specified istream in form parsable by the
       * GeometryInfo::read method.
       */
      void write(std::ostream &oStream) const;
      void writeWithoutVersion(std::ostream &oStream) const;

      bool isCompatible(const GeometryInfo&) const;
      bool isIdenticalGeometry(const GeometryInfo&) const;
      
    private:
      class Impl;
      
      Impl *m_pImpl;
    };
    
    std::ostream &operator<<(std::ostream &oStream, const GeometryInfo &geoInfo);
    std::istream &operator<<(std::istream &iStream, GeometryInfo &geoInfo);
  }
}

#endif
