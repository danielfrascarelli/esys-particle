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



#include <stdexcept>
#include <sstream>
#include <vector>

namespace esys
{
  namespace lsm
  {
    namespace impl
    {
      inline double square(double val)
      {
        return val*val;
      }

      template <int tmplDim, typename TmplVec>
      DimBasicBox<tmplDim,TmplVec>::DimBasicBox(const Vec &minPt, const Vec &maxPt)
        : m_minPt(minPt),
          m_maxPt(maxPt)
      {
      }

      template <int tmplDim, typename TmplVec>
      const typename DimBasicBox<tmplDim,TmplVec>::Vec &DimBasicBox<tmplDim,TmplVec>::getMinPt() const
      {
        return m_minPt;
      }

      template <int tmplDim, typename TmplVec>
      const typename DimBasicBox<tmplDim,TmplVec>::Vec &DimBasicBox<tmplDim,TmplVec>::getMaxPt() const
      {
        return m_maxPt;
      }

      template <int tmplDim, typename TmplVec>
      double DimBasicBox<tmplDim,TmplVec>::getVolume() const
      {
        double product = 1.0;
        for (int i = 0; i < tmplDim; i++)
        {
          product *= (this->getMaxPt()[i] - this->getMinPt()[i]);
        }
        return product;
      }
      
      template <int tmplDim, typename TmplVec>
      template <typename TmplSphere>
      bool DimBasicBox<tmplDim,TmplVec>::intersectsWith(const TmplSphere &sphere) const
      {
        double distSqrd = 0.0;
        for (int i = 0; i < tmplDim; i++)
        {
          if (sphere.getCentre()[i] < this->getMinPt()[i])
          {
            distSqrd += square(sphere.getCentre()[i] - this->getMinPt()[i]);
          }
          else if (sphere.getCentre()[i] > this->getMaxPt()[i])
          {
            distSqrd += square(sphere.getCentre()[i] - this->getMaxPt()[i]);
          }
        }
        return (distSqrd <= square(sphere.getRadius()));
      }

      template <int tmplDim, typename TmplVec>
      bool DimBasicBox<tmplDim,TmplVec>::intersectsWith(const Vec &pt) const
      {
        for (int i = 0; (i < tmplDim); i++)
        {
          if ((this->getMinPt()[i] > pt[i]) || (pt[i] > this->getMaxPt()[i]))
          {
            return false;
          }
        }
        return true;
      }

      template <int tmplDim, typename TmplVec>
      template <typename TmplSphere>
      bool DimBasicBox<tmplDim,TmplVec>::contains(const TmplSphere &sphere) const
      {
        for (int i = 0; (i < tmplDim); i++)
        {
          Vec pt = Vec(0.0);
          pt[i] = sphere.getRadius();
          if
          (
            !(intersectsWith(sphere.getCentre() + pt))
            ||
            !(intersectsWith(sphere.getCentre() - pt))
          )
          {
            return false;
          }
        }
        return true;
      }

      template <int tmplDim, typename TmplVec>
      double DimPlane<tmplDim,TmplVec>::norm(const Vec &pt)
      {
        double sum = 0.0;
        for (int i = 0; i < tmplDim; i++)
        {
          sum += pt[i]*pt[i];
        }
        return sqrt(sum);
      }

      template <int tmplDim, typename TmplVec>
      double DimPlane<tmplDim,TmplVec>::dot(const Vec &p1, const Vec &p2)
      {
        double sum = 0.0;
        for (int i = 0; i < tmplDim; i++)
        {
          sum += p1[i]*p2[i];
        }
        return sum;
      }
      
      template <int tmplDim, typename TmplVec>
      DimPlane<tmplDim,TmplVec>::DimPlane() : m_normal(), m_pt(), m_invNormalNorm(0.0)
      {
      }

      template <int tmplDim, typename TmplVec>
      DimPlane<tmplDim,TmplVec>::DimPlane(const Vec &normal, const Vec &pt)
        : m_normal(normal),
          m_pt(pt),
          m_invNormalNorm((1.0/norm(normal)))
      {
      }

      template <int tmplDim, typename TmplVec>
      DimPlane<tmplDim,TmplVec>::DimPlane(const DimPlane &plane)
        : m_normal(plane.m_normal),
          m_pt(plane.m_pt),
          m_invNormalNorm(plane.m_invNormalNorm)
      {
      }

      template <int tmplDim, typename TmplVec>
      DimPlane<tmplDim,TmplVec> &DimPlane<tmplDim,TmplVec>::operator=(const DimPlane &plane)
      {
        m_normal        = plane.m_normal;
        m_pt            = plane.m_pt;
        m_invNormalNorm = plane.m_invNormalNorm;
        return *this;
      }

      template <int tmplDim, typename TmplVec>
      double DimPlane<tmplDim,TmplVec>::getSignedDistanceTo(const Vec &pt) const
      {
        // http://mathworld.wolfram.com/Point-PlaneDistance.html
        return
          (
            (dot(m_normal, pt) - dot(m_normal, m_pt))*m_invNormalNorm
          );
      }

      template <int tmplDim, typename TmplVec>
      double DimPlane<tmplDim,TmplVec>::getDistanceTo(const Vec &pt) const
      {
        return fabs(getSignedDistanceTo(pt));
      }

      template <int tmplDim, typename TmplVec>
      const typename DimPlane<tmplDim,TmplVec>::Vec &DimPlane<tmplDim,TmplVec>::getNormal() const
      {
        return m_normal;
      }

      template <int tmplDim, typename TmplVec>
      const double DimBasicSphere<tmplDim,TmplVec>::FOUR_THIRDS_PI = (4.0/3.0)*M_PI;

      template <int tmplDim, typename TmplVec>
      const double DimBasicSphere<tmplDim,TmplVec>::ONE_THIRD_PI   = M_PI/3.0;

      template <int tmplDim, typename TmplVec>
      DimBasicSphere<tmplDim,TmplVec>::DimBasicSphere()
        : m_centre(),
          m_radius(0.0)
      {
      }

      template <int tmplDim, typename TmplVec>
      DimBasicSphere<tmplDim,TmplVec>::DimBasicSphere(const Vec &centrePt, double radius)
        : m_centre(centrePt),
          m_radius(radius)
      {
      }

      template <int tmplDim, typename TmplVec>
      DimBasicSphere<tmplDim,TmplVec>::DimBasicSphere(const DimBasicSphere &sphere)
        : m_centre(sphere.m_centre),
          m_radius(sphere.m_radius)
      {
      }

      template <int tmplDim, typename TmplVec>
      DimBasicSphere<tmplDim,TmplVec> &DimBasicSphere<tmplDim,TmplVec>::operator=(const DimBasicSphere &sphere)
      {
        m_centre = sphere.m_centre;
        m_radius = sphere.m_radius;
        return *this;
      }

      template <int tmplDim, typename TmplVec>
      double DimBasicSphere<tmplDim,TmplVec>::getRadius() const
      {
        return m_radius;
      }

      template <int tmplDim, typename TmplVec>
      const typename DimBasicSphere<tmplDim,TmplVec>::Vec &
      DimBasicSphere<tmplDim,TmplVec>::getCentre() const
      {
        return m_centre;
      }

      template <int tmplDim, typename TmplVec>
      double DimBasicSphere<tmplDim,TmplVec>::getVolume() const
      {
        return (tmplDim == 2) ? M_PI*getRadius()*getRadius() : FOUR_THIRDS_PI*getRadius()*getRadius()*getRadius();
      }

      inline void checkDomain(double r, double x1, double y1, double x2, double y2)
      {
        const double rSqrd = r*r;
        const double x1Sqrd = x1*x1;
        const double x2Sqrd = x2*x2;
        const double y1Sqrd = y1*y1;
        const double y2Sqrd = y2*y2;
        if
        (
          ((rSqrd - x1Sqrd - y1Sqrd) < 0)
          ||
          ((rSqrd - x1Sqrd - y2Sqrd) < 0)
          ||
          ((rSqrd - x2Sqrd - y1Sqrd) < 0)
          ||
          ((rSqrd - x2Sqrd - y2Sqrd) < 0)
        )
        {
          std::stringstream msg;
          msg 
            << "Invalid rectangular domain for sphere integration, points in domain "
            << "(" << x1 << "," << y1 << "), (" << x2 << "," << y2 << " lie outside "
            << "sphere of radius " << r << " centred at the origin.";
          throw std::runtime_error(msg.str());
        }
      }

      template <int tmplDim, typename TmplVec>
      double DimBasicSphere<tmplDim,TmplVec>::getVolume(const Vec &minPt, const Vec &maxPt, const int dimX, const int dimY) const
      {
        double vol = 0.0;
        if ((tmplDim == 2) || (tmplDim == 3))
        {
          if (minPt[dimX] != maxPt[dimX])
          {
            const double x1 = minPt[dimX] - getCentre()[dimX];
            const double x2 = maxPt[dimX] - getCentre()[dimX];
            const double r  = getRadius();
  
            if (tmplDim == 2)
            {
              const double rSqrd  = r*r;
              const double x1Sqrd = x1*x1;
              const double x2Sqrd = x2*x2;

              vol = 
                (
                  rSqrd*asin(x2/r)
                  +
                  x2*sqrt(rSqrd-x2Sqrd)
                  -
                  rSqrd*asin(x1/r)
                  -
                  x1*sqrt(rSqrd-x1Sqrd)
                )*0.5;
            }
            else if (tmplDim == 3)
            {
              if (minPt[dimY] != maxPt[dimY])
              {
                const double y1 = minPt[dimY] - getCentre()[dimY];
                const double y2 = maxPt[dimY] - getCentre()[dimY];

                checkDomain(r, x1, y1, x2, y2);

                //
                // Matlab6/maple generated code:
                //
                // syms x y x1 x2 y1 y2 r real
                // sphereIntegral = int(int('sqrt(r^2-x^2-y^2)', x, x1, x2), y, y1, y2)
                // maple('readlib(codegen)')
                // maple('readlib(optimize)')
                // optimizedIntegral = maple('optimize', sphereIntegral, 'tryhard')
                // maple('cost', optimizedIntegral)
                //
                //43*additions+82*multiplications+6*divisions+22*functions+42*assignments

                const double t30 = y2*y2; //y2^2;
                const double t31 = x2*x2; //x2^2;
                const double t36 = r*r;   //r^2;
                const double t59 = t31-t36;
                const double t40 = sqrt(-t30-t59); //pow(-t30-t59,1/2);
                const double t10 = 1.0/t40;
                const double t32 = x1*x1; //x1^2;
                const double t54 = t32-t36;
                const double t42 = sqrt(-t30-t54); //pow(-t30-t54,1/2);
                const double t14 = 1.0/t42;
                const double t64 = -atan(t10*x2)+atan(t14*x1);
                const double t27 = y1*y1; //y1^2;
                const double t39 = sqrt(-t27-t59); //pow(-t27-t59,1/2);
                const double t9  = 1.0/t39;
                const double t41 = sqrt(-t27-t54); //pow(-t27-t54,1/2);
                const double t12 = 1.0/t41;
                const double t63 = atan(t12*x1)-atan(t9*x2);
                const double t62 = -atan(y2*t14)+atan(y1*t12);
                const double t61 = -atan(y1*t9)+atan(t10*y2);
                const double t37 = sqrt(t31); //pow(t31,1/2);
                const double t21 = 1.0/t37;
                const double t60 = t21*t9;
                const double t38 = sqrt(t32); //pow(t32,1/2);
                const double t24 = 1.0/t38;
                const double t58 = t24*t14;
                const double t57 = t24*t12;
                const double t56 = t37*t38;
                const double t55 = t21*t10;
                const double t53 = 2.0*x2;
                const double t52 = 2.0*x1;
                const double t51 = t42*t56;
                const double t28 = t27*y1;
                const double t50 = t28-t36*y1;
                const double t34 = t30*y2;
                const double t49 = t34-t36*y2;
                const double t48 = t41*t51;
                const double t35 = t36*r;
                const double t33 = t31*x2;
                const double t29 = t32*x1;
                const double t26 = r*y2;
                const double t25 = y1*r;
                vol = 
                  (-1.0/6.0)
                  *
                  (
                    (-2.0*t33*y1-t50*t53)*t40*t48
                    +
                    (
                      (2.0*t33*y2+t49*t53)*t48
                      +
                      (
                        (2.0*t29*y1+t50*t52)*t51
                        +
                        (
                          (-2.0*t29*y2-t49*t52)*t56
                          +
                          (
                            (
                              atan((t25+t54)*t57)
                              +
                              atan((-t26+t54)*t58)
                              -
                              atan((t26+t54)*t58)
                              -
                              atan((-t25+t54)*t57)
                            )*t35*t37*x1
                            +
                            (
                              (
                                -atan((-t26+t59)*t55)
                                -
                                atan((t25+t59)*t60)
                                +
                                atan((-t25+t59)*t60)
                                +
                                atan((t26+t59)*t55)
                              )*t35*x2
                              +
                              (-t64*t34+t61*t33+t62*t29+t63*t28+3.0*(t64*y2-t63*y1-t61*x2-t62*x1)*t36)*t37
                            )*t38
                          )*t42
                        )*t41
                      )*t40
                    )*t39
                  )*t14*t9*t55*t57;
              }
            }
          }
        }
        
        return vol;
      }

      template <int tmplDim, typename TmplVec>
      bool DimBasicSphere<tmplDim,TmplVec>::intersectsWith(const Vec &pt) const
      {
        double distSqrd = 0.0;
        for (int i = 0; i < tmplDim; i++)
        {
          distSqrd += square(getCentre()[i] - pt[i]);
        }
        return (distSqrd <= square(getRadius()));
      }

      template <int tmplDim, typename TmplVec>
      double DimBasicSphere<tmplDim,TmplVec>::getSegmentVolume(const Plane &plane) const
      {
        double vol = 0.0;
        if ((tmplDim == 2) || (tmplDim == 3))
        {
          const double signedD = plane.getSignedDistanceTo(getCentre());
          const double d = fabs(signedD);
          if (d < getRadius())
          {
            if (tmplDim == 2)
            {
              // http://mathworld.wolfram.com/CircularSegment.html
              const double rSqrd = getRadius()*getRadius();
              vol = rSqrd*acos(d/getRadius()) - d*sqrt(rSqrd - d*d);
            }
            else if (tmplDim == 3)
            {
              // http://mathworld.wolfram.com/SphericalCap.html
              const double h = getRadius() - d;
              vol = ONE_THIRD_PI*h*h*(3.0*getRadius()-h);
            }
            vol = ((signedD < 0) ? vol : getVolume() - vol);
          }
        }
        return vol;
      }

      template <int tmplDim, typename TmplVec>
      typename IntersectionVolCalculator<tmplDim,TmplVec>::Vec
      IntersectionVolCalculator<tmplDim,TmplVec>::getNormal(int dim)
      {
        Vec n = Vec(0.0);
        n[dim] = 1.0;
        return n;
      }

      template <int tmplDim, typename TmplVec>
      typename IntersectionVolCalculator<tmplDim,TmplVec>::Vec
      IntersectionVolCalculator<tmplDim,TmplVec>::getNegNormal(int dim)
      {
        Vec n = Vec(0.0);
        n[dim] = -1.0;
        return n;
      }

      template <int tmplDim, typename TmplVec>
      IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere::VolumeSphere()
        : m_sphere(),
          m_volume(0.0)
      {
      }

      template <int tmplDim, typename TmplVec>
      IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere::VolumeSphere(
        const BasicSphere &sphere
      )
        : m_sphere(sphere),
          m_volume(sphere.getVolume())
      {
      }

      template <int tmplDim, typename TmplVec>
      IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere::VolumeSphere(
        const VolumeSphere &sphere
      )
        : m_sphere(BasicSphere(sphere.getCentre(), sphere.getRadius())),
          m_volume(sphere.m_volume)
      {
      }

      template <int tmplDim, typename TmplVec>
      typename IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere &
      IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere::operator=(
        const VolumeSphere &sphere
      )
      {
        m_sphere = sphere.m_sphere;
        m_volume = sphere.m_volume;
        return *this;
      }

      template <int tmplDim, typename TmplVec>
      double IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere::getRadius() const
      {
        return m_sphere.getRadius();
      }

      template <int tmplDim, typename TmplVec>
      const typename IntersectionVolCalculator<tmplDim,TmplVec>::Vec &
      IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere::getCentre() const
      {
        return m_sphere.getCentre();
      }

      template <int tmplDim, typename TmplVec>
      double IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere::getVolume() const
      {
        return m_volume;
      }

      template <int tmplDim, typename TmplVec>
      double IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere::getVolume(
        const Vec &minPt,
        const Vec &maxPt,
        const int dimX,
        const int dimY
      ) const
      {
        return m_sphere.getVolume(minPt, maxPt, dimX, dimY);
      }

      template <int tmplDim, typename TmplVec>
      double IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere::calcVolume() const
      {
        return m_sphere.getVolume();
      }

      template <int tmplDim, typename TmplVec>
      bool IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere::intersectsWith(const Vec &pt) const
      {
        return m_sphere.intersectsWith(pt);
      }

      template <int tmplDim, typename TmplVec>
      double IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere::getSegmentVolume(const Plane &plane) const
      {
        return m_sphere.getSegmentVolume(plane);
      }

      template <int tmplDim, typename TmplVec>
      IntersectionVolCalculator<tmplDim,TmplVec>::Vertex::Vertex() : m_pt()
      {
      }

      template <int tmplDim, typename TmplVec>
      IntersectionVolCalculator<tmplDim,TmplVec>::Vertex::Vertex(const Vec &pt) : m_pt(pt)
      {
      }

      template <int tmplDim, typename TmplVec>
      IntersectionVolCalculator<tmplDim,TmplVec>::Vertex::Vertex(const Vertex &vtx) : m_pt(vtx.m_pt)
      {
      }

      template <int tmplDim, typename TmplVec>
      typename IntersectionVolCalculator<tmplDim,TmplVec>::Vertex &
      IntersectionVolCalculator<tmplDim,TmplVec>::Vertex::operator=(const Vertex &vtx)
      {
        m_pt = vtx.m_pt;
        return *this;
      }

      template <int tmplDim, typename TmplVec>
      const typename IntersectionVolCalculator<tmplDim,TmplVec>::Vec &
      IntersectionVolCalculator<tmplDim,TmplVec>::Vertex::getPoint() const
      {
        return m_pt;
      }

      template <int tmplDim, typename TmplVec>
      void IntersectionVolCalculator<tmplDim,TmplVec>::Vertex::setPoint(const Vec &pt)
      {
        m_pt = pt;
      }

      template <int tmplDim, typename TmplVec>
      IntersectionVolCalculator<tmplDim,TmplVec>::VertexBox::VertexBox(
        const BasicBox &box
      )
        : BasicBox(box),
          m_vertexArray()
      {
        createVertices();
      }

      template <int tmplDim, typename TmplVec>
      IntersectionVolCalculator<tmplDim,TmplVec>::VertexBox::VertexBox(
        const VertexBox &box
      )
        : BasicBox(box)
      {
        for (int i = 0; i < getNumVertices(); i++)
        {
          m_vertexArray[i] = box.m_vertexArray[i];
        }
      }

      template <int tmplDim, typename TmplVec>
      typename IntersectionVolCalculator<tmplDim,TmplVec>::VertexBox &
      IntersectionVolCalculator<tmplDim,TmplVec>::VertexBox::operator=(
        const VertexBox &box
      )
      {
        BasicBox::operator=(box);
        for (int i = 0; i < getNumVertices(); i++)
        {
          m_vertexArray[i] = box.m_vertexArray[i];
        }
        return *this;
      }

      template <int tmplDim, typename TmplVec>
      void IntersectionVolCalculator<tmplDim,TmplVec>::VertexBox::createVertices()
      {
        int j = 0;
        m_vertexArray[j].setPoint(this->getMinPt());
        int i = 0;
        for (j++; i < tmplDim; i++, j++)
        {
          Vec pt = this->getMinPt();
          pt[i] = this->getMaxPt()[i];
          m_vertexArray[j].setPoint(pt);
        }

        m_vertexArray[j].setPoint(this->getMaxPt());
        for (i = 0, j++; i < tmplDim && j < s_numVertices; i++, j++)
        {
          Vec pt = this->getMaxPt();
          pt[i] = this->getMinPt()[i];
          m_vertexArray[j] = pt;
        }
      }

      template <int tmplDim, typename TmplVec>
      const typename IntersectionVolCalculator<tmplDim,TmplVec>::Vertex &
      IntersectionVolCalculator<tmplDim,TmplVec>::VertexBox::getVertex(int i) const
      {
        return m_vertexArray[i];
      }

      template <int tmplDim, typename TmplVec>
      int IntersectionVolCalculator<tmplDim,TmplVec>::VertexBox::getNumVertices()
      {
        return s_numVertices;
      }

      template <int tmplDim, typename TmplVec>
      IntersectionVolCalculator<tmplDim,TmplVec>::IntersectionVolCalculator(
        const BasicBox &box
      )
        : m_sphere(BasicSphere(Vec(), 1.0)),
          m_box(box)
      {
      }

      template <int tmplDim, typename TmplVec>
      const typename IntersectionVolCalculator<tmplDim,TmplVec>::VolumeSphere &
      IntersectionVolCalculator<tmplDim,TmplVec>::getSphere() const
      {
        return m_sphere;
      }

      template <int tmplDim, typename TmplVec>
      void IntersectionVolCalculator<tmplDim,TmplVec>::setSphere(
        const BasicSphere &sphere
      )
      {
        m_sphere = sphere;
      }

      template <int tmplDim, typename TmplVec>
      const typename IntersectionVolCalculator<tmplDim,TmplVec>::BasicBox &
      IntersectionVolCalculator<tmplDim,TmplVec>::getBox() const
      {
        return m_box;
      }

      template <int tmplDim, typename TmplVec>
      const typename IntersectionVolCalculator<tmplDim,TmplVec>::VertexBox &
      IntersectionVolCalculator<tmplDim,TmplVec>::getVertexBox() const
      {
        return m_box;
      }

      template <int tmplDim, typename TmplVec>
      typename IntersectionVolCalculator<tmplDim,TmplVec>::Vec
      IntersectionVolCalculator<tmplDim,TmplVec>::componentMin(
        const Vec &p1,
        const Vec &p2
      )
      {
        Vec m;
        for (int i = 0; i < tmplDim; i++)
        {
          m[i] = min(p1[i], p2[i]);
        }
        return m;
      }

      template <int tmplDim, typename TmplVec>
      typename IntersectionVolCalculator<tmplDim,TmplVec>::Vec
      IntersectionVolCalculator<tmplDim,TmplVec>::componentMax(
        const Vec &p1,
        const Vec &p2
      )
      {
        Vec m;
        for (int i = 0; i < tmplDim; i++)
        {
          m[i] = max(p1[i], p2[i]);
        }
        return m;
      }

      template <int tmplDim, typename TmplVec>
      double IntersectionVolCalculator<tmplDim,TmplVec>::getInsidePointVolume(
        const Vec &pt
      ) const
      {
        double vol = 0.0;
        const Vec &centrePt = getSphere().getCentre();
        const Vec diag = (getSphere().getCentre() - pt)*2.0;
        const Vec oppCorner = pt + diag;
        BasicBox box =
          BasicBox(
            componentMin(pt, oppCorner),
            componentMax(pt, oppCorner)
          );
        const double boxVol  = box.getVolume();
        const double sphVol  = getSphere().getVolume();

        double s[tmplDim];
        double v[tmplDim+1];
        for (int i = 0; i < tmplDim; i++)
        {
          s[i] = getSphere().getSegmentVolume(Plane(getNormal((i+1)%tmplDim), box.getMaxPt()));
        }
        if (tmplDim == 2)
        {
          v[0] = 0.50*(sphVol - 2.0*s[0] - boxVol);
          v[1] = 0.50*(sphVol - 2.0*s[1] - boxVol);
          v[2] = 0.25*(sphVol - 2.0*v[0] - 2.0*v[1] - boxVol);

          if (pt[0] <= centrePt[0])
          {
            if (pt[1] <= centrePt[1])
            {
              vol = boxVol + v[0] + v[1] + v[2];
            }
            else
            {
              vol = v[1] + v[2];
            }
          }
          else
          {
            if (pt[1] <= centrePt[1])
            {
              vol = v[0] + v[2];
            }
            else
            {
              vol = v[2];
            }
          }
        }
        else if (tmplDim == 3)
        {
          v[0] = 
            0.500*(
              2.0*getSphere().getVolume(
                box.getMinPt(),
                box.getMaxPt(),
                1,
                2
              )
              -
              boxVol
            );
          v[1] = 
            0.500*(
              2.0*getSphere().getVolume(
                box.getMinPt(),
                box.getMaxPt(),
                0,
                2
              )
              -
              boxVol
            );
          v[2] = 
            0.500*(
              2.0*getSphere().getVolume(
                box.getMinPt(),
                box.getMaxPt(),
                0,
                1
              )
              -
              boxVol
            );

          double e[3];
          e[0] = 0.250*(sphVol - 2.0*s[1] - 2.0*v[0] - 2.0*v[1] - boxVol);
          e[1] = 0.250*(sphVol - 2.0*s[2] - 2.0*v[1] - 2.0*v[2] - boxVol);
          e[2] = 0.250*(sphVol - 2.0*s[0] - 2.0*v[0] - 2.0*v[2] - boxVol);
          
          v[3] = 
            0.125*(
              sphVol
              -
              2.0*v[0] - 2.0*v[1] - 2.0*v[2]
              -
              4.0*e[0] - 4.0*e[1] - 4.0*e[2]
              -
              boxVol
            );

          if (pt[0] <= centrePt[0])
          {
            if (pt[1] <= centrePt[1])
            {
              if (pt[2] <= centrePt[2])
              {
                vol = boxVol + v[0] + v[1] + v[2] + v[3] + e[0] + e[1] + e[2];
              }
              else
              {
                vol = v[2] + v[3] + e[1] + e[2];
              }
            }
            else
            {
              if (pt[2] <= centrePt[2])
              {
                vol = v[1] + v[3] + e[0] + e[1];
              }
              else
              {
                vol = v[3] + e[1];
              }
            }
          }
          else
          {
            if (pt[1] <= centrePt[1])
            {
              if (pt[2] <= centrePt[2])
              {
                vol = v[0] + v[3] + e[0] + e[2];
              }
              else
              {
                vol = v[3] + e[2];
              }
            }
            else
            {
              if (pt[2] <= centrePt[2])
              {
                vol = v[3] + e[0];
              }
              else
              {
                vol = v[3];
              }
            }
          }
        }

        return vol;
      }

      template <int tmplDim, typename TmplVec>
      double IntersectionVolCalculator<tmplDim,TmplVec>::getTwoPlaneVolume(
        const Vec &pt,
        const int orientDim
      ) const
      {
        const int ZERO = orientDim;
        const int ONE  = (orientDim+1) % tmplDim;
        const int TWO  = (orientDim+2) % tmplDim;

        double vol = 0.0;
        const double sphVol = getSphere().getVolume();
        const Vec &centrePt = getSphere().getCentre();

        if ((square(pt[ONE]-centrePt[ONE]) + square(pt[TWO]-centrePt[TWO])) < square(getSphere().getRadius()))
        {
          Plane plane[tmplDim];
          plane[ZERO] = Plane();
          plane[ONE]  = Plane(getNormal(ONE), pt);
          plane[TWO]  = Plane(getNormal(TWO), pt);

          const double halfSphVol = sphVol*0.5;
          double s[tmplDim];
          s[ZERO] = 0;
          s[ONE]  = getSphere().getSegmentVolume(plane[ONE]);
          s[TWO]  = getSphere().getSegmentVolume(plane[TWO]);
          s[ONE] = ((s[ONE] > halfSphVol) ? (sphVol - s[ONE]) : s[ONE]);
          s[TWO] = ((s[TWO] > halfSphVol) ? (sphVol - s[TWO]) : s[TWO]);

          Vec distVec(4.0*getSphere().getRadius());
          distVec[ONE] = plane[ONE].getDistanceTo(centrePt);
          distVec[TWO] = plane[TWO].getDistanceTo(centrePt);
          const double coreVol =
            2.0*getSphere().getVolume(
              centrePt - Vec(distVec[ONE], distVec[TWO], distVec[ZERO]),
              centrePt + Vec(distVec[ONE], distVec[TWO], distVec[ZERO])
            );
          double v[tmplDim];
          v[ONE]  = 0.50*(sphVol - 2.0*s[TWO] - coreVol);
          v[TWO]  = 0.50*(sphVol - 2.0*s[ONE] - coreVol);
          v[ZERO] = 0.25*(sphVol - 2.0*v[ONE] - 2.0*v[TWO] - coreVol);

          if (pt[ONE] <= centrePt[ONE])
          {
            if (pt[TWO] <= centrePt[TWO])
            {
              vol = coreVol + v[ONE] + v[TWO] + v[ZERO];
            }
            else
            {
              vol = v[TWO] + v[ZERO];
            }
          }
          else
          {
            if (pt[TWO] <= centrePt[TWO])
            {
              vol = v[ONE] + v[ZERO];
            }
            else
            {
              vol = v[ZERO];
            }
          }
        }
        else
        {
          if (pt[ONE] <= centrePt[ONE])
          {
            if (pt[TWO] <= centrePt[TWO])
            {
              vol = 
                sphVol
                -
                getSphere().getSegmentVolume(Plane(getNegNormal(ONE), pt))
                -
                getSphere().getSegmentVolume(Plane(getNegNormal(TWO), pt));
            }
            else
            {
              vol = getSphere().getSegmentVolume(Plane(getNormal(TWO), pt));
            }
          }
          else
          {
            if (pt[TWO] <= centrePt[TWO])
            {
              vol = getSphere().getSegmentVolume(Plane(getNormal(ONE), pt));
            }
            else
            {
              vol = 0.0;
            }
          }
        }

        return vol;
      }

      template <int tmplDim, typename TmplVec>
      double IntersectionVolCalculator<tmplDim,TmplVec>::getOutsidePointVolume(
        const Vec &pt
      ) const
      {
        double vol = 0.0;
        const double sphVol = getSphere().getVolume();
        const Vec &centrePt = getSphere().getCentre();
        if (tmplDim == 2)
        {
          if (pt[0] <= centrePt[0])
          {
            if (pt[1] <= centrePt[1])
            {
              vol = 
                sphVol
                -
                getSphere().getSegmentVolume(Plane(getNegNormal(0), pt))
                -
                getSphere().getSegmentVolume(Plane(getNegNormal(1), pt));
            }
            else
            {
              vol = getSphere().getSegmentVolume(Plane(getNormal(1), pt));
            }
          }
          else
          {
            if (pt[1] <= centrePt[1])
            {
              vol = getSphere().getSegmentVolume(Plane(getNormal(0), pt));
            }
            else
            {
              vol = 0.0;
            }
          }
        }
        else if (tmplDim == 3)
        {
          const Vec diag = (centrePt - pt)*2.0;
          const Vec oppCorner = pt + diag;
          BasicBox box =
            BasicBox(
              componentMin(pt, oppCorner),
              componentMax(pt, oppCorner)
            );

          double s[tmplDim];
          double e[tmplDim];
          for (int i = 0; i < tmplDim; i++)
          {
            s[i] = getSphere().getSegmentVolume(Plane(getNormal(i), box.getMaxPt()));
            e[i] = getTwoPlaneVolume(box.getMaxPt(), i);
          }
          double v[tmplDim+1];
          v[0] = s[0] - 2.0*e[1] - 2.0*e[2];
          v[1] = s[1] - 2.0*e[0] - 2.0*e[2];
          v[2] = s[2] - 2.0*e[0] - 2.0*e[1];
          v[3] = sphVol - (4.0*e[0] + 4.0*e[1] + 4.0*e[2] + 2.0*v[0] + 2.0*v[1] + 2.0*v[2]);

          if (pt[0] <= centrePt[0])
          {
            if (pt[1] <= centrePt[1])
            {
              if (pt[2] <= centrePt[2])
              {
                vol = v[0] + v[1] + v[2] + v[3] + e[0] + e[1] + e[2];
              }
              else
              {
                vol = v[2] + e[0] + e[1];
              }
            }
            else
            {
              if (pt[2] <= centrePt[2])
              {
                vol = v[1] + e[0] + e[2];
              }
              else
              {
                vol = e[0];
              }
            }
          }
          else
          {
            if (pt[1] <= centrePt[1])
            {
              if (pt[2] <= centrePt[2])
              {
                vol = v[0] + e[1] + e[2];
              }
              else
              {
                vol = e[1];
              }
            }
            else
            {
              if (pt[2] <= centrePt[2])
              {
                vol = e[2];
              }
              else
              {
                vol = 0.0;
              }
            }
          }
        }
        return vol;
      }

      
      template <int tmplDim, typename TmplVec>
      double IntersectionVolCalculator<tmplDim,TmplVec>::getVolume(
        const Vertex &vtx
      )
      {
        double vol = 0;
        if (getSphere().intersectsWith(vtx.getPoint()))
        {
          vol = getInsidePointVolume(vtx.getPoint());
        }
        else
        {
          vol = getOutsidePointVolume(vtx.getPoint());
        }
        return vol;
      }

      template <int tmplDim, typename TmplVec>
      double IntersectionVolCalculator<tmplDim,TmplVec>::getVertexVolume(
        const BasicSphere &sphere
      )
      {
        m_sphere = VolumeSphere(sphere);
        double vol = 0.0;

	std::vector<double> vtxVol(getVertexBox().getNumVertices());
        for (int i = 0; i < getVertexBox().getNumVertices(); i++)
        {
          vtxVol[i] = getVolume(getVertexBox().getVertex(i));
        }

        if (tmplDim == 2)
        {
          vol = vtxVol[0] - vtxVol[1] - vtxVol[2] + vtxVol[3];
        }
        else if (tmplDim == 3)
        {
          vol = 
            vtxVol[7] + vtxVol[6] + vtxVol[5] - vtxVol[4]
            -
            vtxVol[3] - vtxVol[2] - vtxVol[1] + vtxVol[0];
        }

        return vol;
      }

      template <int tmplDim, typename TmplVec>
      bool IntersectionVolCalculator<tmplDim,TmplVec>::sphereContainsBox(
        const BasicSphere &sphere
      ) const
      {
        for (int i = 0; i < getVertexBox().getNumVertices(); i++)
        {
          if (!sphere.intersectsWith(getVertexBox().getVertex(i).getPoint()))
          {
            return false;
          }
        }
        return true;
      }

      template <int tmplDim, typename TmplVec>
      double IntersectionVolCalculator<tmplDim,TmplVec>::getVolume(
        const BasicSphere &sphere
      )
      {
        double vol = 0.0;
        if (getBox().intersectsWith(sphere))
        {
          if (sphereContainsBox(sphere))
          {
            vol = getBox().getVolume();
          }
          else if (getBox().contains(sphere))
          {
            vol = sphere.getVolume();
          }
          else
          {
            vol = getVertexVolume(sphere);
          }
        }
        return vol;
      }
    }
  }
}
