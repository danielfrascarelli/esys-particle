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


#ifndef ESYS_LSM_IMPLINTERSECTIONVOLCALCULATOR_H
#define ESYS_LSM_IMPLINTERSECTIONVOLCALCULATOR_H

#include <math.h>

namespace esys
{
  namespace lsm
  {
    namespace impl
    {
      double square(double val);

      template <int tmplDim, typename TmplVec>
      class DimBasicBox
      {
      public:
        typedef TmplVec Vec;
        DimBasicBox(const Vec &minPt, const Vec &maxPt);
  
        const Vec &getMinPt() const;
  
        const Vec &getMaxPt() const;
  
        double getVolume() const;
        
        template <typename TmplSphere>
        bool intersectsWith(const TmplSphere &sphere) const;
  
        bool intersectsWith(const Vec &pt) const;
  
        template <typename TmplSphere>
        bool contains(const TmplSphere &sphere) const;
  
      private:
        Vec m_minPt;
        Vec m_maxPt;
      };
  
      template <int tmplDim, typename TmplVec>
      class DimPlane
      {
      public:
        typedef TmplVec Vec;
  
        static double norm(const Vec &pt);
  
        static double dot(const Vec &p1, const Vec &p2);
        
        DimPlane();
  
        DimPlane(const Vec &normal, const Vec &pt);
  
        DimPlane(const DimPlane &plane);
  
        DimPlane &operator=(const DimPlane &plane);
  
        double getSignedDistanceTo(const Vec &pt) const;
  
        double getDistanceTo(const Vec &pt) const;
  
        const Vec &getNormal() const;
  
      private:
        Vec    m_normal;
        Vec    m_pt;
        double m_invNormalNorm;
      };
  
      template <int tmplDim, typename TmplVec>
      class DimBasicSphere
      {
      public:
        typedef TmplVec                 Vec;
        typedef DimPlane<tmplDim, Vec> Plane;
  
        static const double FOUR_THIRDS_PI;
        static const double ONE_THIRD_PI;

        DimBasicSphere();
        
        DimBasicSphere(const Vec &centrePt, double radius);

        DimBasicSphere(const DimBasicSphere &sphere);

        DimBasicSphere &operator=(const DimBasicSphere &sphere);

        double getRadius() const;

        const Vec &getCentre() const;

        double getVolume() const;

        double getVolume(const Vec &minPt, const Vec &maxPt, const int dimX = 0, const int dimY = 1) const;

        bool intersectsWith(const Vec &pt) const;

        double getSegmentVolume(const Plane &plane) const;

      private:
        Vec    m_centre;
        double m_radius;
      };
  
      template <int tmplDim, typename TmplVec>
      class IntersectionVolCalculator
      {
      public:
        typedef TmplVec                     Vec;
        typedef DimBasicSphere<tmplDim,Vec> BasicSphere;
        typedef DimBasicBox<tmplDim,Vec>    BasicBox;
        typedef DimPlane<tmplDim,Vec>       Plane;
  
        static Vec getNormal(int dim);
  
        static Vec getNegNormal(int dim);
        
        class VolumeSphere
        {
        public:
          VolumeSphere();

          VolumeSphere(const BasicSphere &sphere);

          VolumeSphere(const VolumeSphere &sphere);

          VolumeSphere &operator=(const VolumeSphere &sphere);

          double getRadius() const;
          
          const Vec &getCentre() const;
          
          double getVolume() const;

          double getVolume(const Vec &minPt, const Vec &maxPt, const int dimX = 0, const int dimY = 1) const;

          double calcVolume() const;

          bool intersectsWith(const Vec &pt) const;

          double getSegmentVolume(const Plane &plane) const;

        private:
          BasicSphere m_sphere;
          double      m_volume;
        };
  
        class Vertex
        {
        public:
          Vertex();
  
          Vertex(const Vec &pt);
          
          Vertex(const Vertex &vtx);
          
          Vertex &operator=(const Vertex &vtx);
  
          const Vec &getPoint() const;
  
          void setPoint(const Vec &pt);
  
        private:
          Vec        m_pt;
        };
        
        class VertexBox : public BasicBox
        {
        public:
          VertexBox(const BasicBox &box);

          VertexBox(const VertexBox &box);

          VertexBox &operator=(const VertexBox &box);

          void createVertices();

          const Vertex &getVertex(int i) const;

          static int getNumVertices();

        private:
          static const int  s_numVertices = ((tmplDim == 2) ? 4 :  8);
          Vertex            m_vertexArray[s_numVertices];
        };

        IntersectionVolCalculator(const BasicBox &box);

        const VolumeSphere &getSphere() const;

        void setSphere(const BasicSphere &sphere);

        const BasicBox &getBox() const;

        const VertexBox &getVertexBox() const;

        static Vec componentMin(const Vec &p1, const Vec &p2);

        static Vec componentMax(const Vec &p1, const Vec &p2);

        double getInsidePointVolume(const Vec &pt) const;

        double getTwoPlaneVolume(const Vec &pt, const int orientDim) const;

        double getOutsidePointVolume(const Vec &pt) const;

        double getVolume(const Vertex &vtx);

        double getVertexVolume(const BasicSphere &sphere);

        bool sphereContainsBox(const BasicSphere &sphere) const;

        double getVolume(const BasicSphere &sphere);

      private:
        VolumeSphere m_sphere;
        VertexBox    m_box;
      };
    }
  }
}

#include "Geometry/IntersectionVolCalculator.hpp"

#endif
