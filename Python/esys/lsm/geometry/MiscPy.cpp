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


#include <boost/version.hpp>
#include <boost/python.hpp>
#include "Python/esys/lsm/geometry/MiscPy.h"
#include "Python/esys/lsm/geometry/SimpleSpherePy.h"
#include "Python/esys/lsm/util/BoundingBoxPy.h"
#include "Python/esys/lsm/util/Vec3Py.h"

namespace esys
{
  namespace lsm
  {
    class SolidBoxPy : public BoundingBoxPy
    {
    public:
      SolidBoxPy() : BoundingBoxPy()
      {
      }
      
      SolidBoxPy(const Vec3Py &minPt, const Vec3Py &maxPt)
        : BoundingBoxPy(minPt, maxPt)
      {
      }

      SolidBoxPy(const SolidBoxPy &box)
        : BoundingBoxPy(box.getMinPt(), box.getMaxPt())
      {
      }

      SolidBoxPy(const BoundingBoxPy &box)
        : BoundingBoxPy(box.getMinPt(), box.getMaxPt())
      {
      }

      SolidBoxPy(const BoundingBox &box)
        : BoundingBoxPy(box.getMinPt(), box.getMaxPt())
      {
      }

      /**
       * Arvo's algorithm for axis-aligned solid-box/solid-sphere intersection.
       */
      bool intersectsWithSphere(const Vec3 &centre, double radius) const
      {
        double distSqrd = 0.0;
        for (int i = 0; i < 3; i++)
        {
          if (centre[i] < getMinPt()[i])
          {
            const double diff = centre[i] - getMinPt()[i];
            distSqrd += diff*diff;
          }
          else if (centre[i] > getMaxPt()[i])
          {
            const double diff = centre[i] - getMaxPt()[i];
            distSqrd += diff*diff;
          }
        }
        return (distSqrd <= (radius*radius));
      }

      bool intersectsWithSpherePy(const Vec3Py &centre, double radius) const
      {
        return intersectsWithSphere(centre, radius);
      }
            
      bool intersects(const SimpleSpherePy &p) const
      {
        return intersectsWithSphere(p.getPos(), p.getRad());
      }
    };

    class HollowBoxPy : public BoundingBoxPy
    {
    public:
      HollowBoxPy() : BoundingBoxPy()
      {
      }

      HollowBoxPy(const Vec3Py &minPt, const Vec3Py &maxPt)
        : BoundingBoxPy(minPt, maxPt)
      {
      }

      HollowBoxPy(const HollowBoxPy &box)
        : BoundingBoxPy(box.getMinPt(), box.getMaxPt())
      {
      }

      HollowBoxPy(const BoundingBoxPy &box)
        : BoundingBoxPy(box.getMinPt(), box.getMaxPt())
      {
      }

      HollowBoxPy(const BoundingBox &box)
        : BoundingBoxPy(box.getMinPt(), box.getMaxPt())
      {
      }

      /**
       * Arvo's algorithm for axis-aligned hollow-box/solid-sphere intersection.
       */
      bool intersectsWithSphere(const Vec3 &centre, double radius) const
      {
        double distSqrd = 0;
        bool face = false;
        for(int i = 0; i < 3; i++)
        {
          if(centre[i] < getMinPt()[i])
          {
            face = true;
            const double diff = (centre[i] - getMinPt()[i]);
            distSqrd += diff*diff;
          }
          else if (centre[i] > getMaxPt()[i])
          {
            face = true;
            const double diff = centre[i] - getMaxPt()[i];
            distSqrd += diff*diff;
          }
          else if ((centre[i] - getMinPt()[i]) <= radius )
          {
            face = true;
          }
          else if ((getMaxPt()[i] - centre[i]) <= radius)
          {
            face = true;
          }
        }
        return (face && (distSqrd <= radius*radius));
      }

      bool intersectsWithSpherePy(const Vec3Py &centre, double radius) const
      {
        return intersectsWithSphere(centre, radius);
      }
      
      bool intersects(const SimpleSpherePy &p) const
      {
        return intersectsWithSphere(p.getPos(), p.getRad());
      }
    };
    
    void exportMisc()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<SolidBoxPy, boost::python::bases<BoundingBoxPy> >(
        "SolidBox"
      )
      .def(boost::python::init<>())
      .def(boost::python::init<const Vec3Py &, const Vec3Py &>())
      .def(boost::python::init<const SolidBoxPy &>())
      .def(boost::python::init<const BoundingBoxPy &>())
      .def("intersects", &SolidBoxPy::intersects)
      .def("intersects", &SolidBoxPy::intersectsWithSpherePy)
      ;

      boost::python::class_<HollowBoxPy, boost::python::bases<BoundingBoxPy> >(
        "HollowBox"
      )
      .def(boost::python::init<>())
      .def(boost::python::init<const Vec3Py &, const Vec3Py &>())
      .def(boost::python::init<const HollowBoxPy &>())
      .def(boost::python::init<const BoundingBoxPy &>())
      .def("intersects", &HollowBoxPy::intersects)
      .def("intersects", &HollowBoxPy::intersectsWithSpherePy)
      ;
    }
  }
}
