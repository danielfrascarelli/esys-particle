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
#include <boost/noncopyable.hpp>
#include "Python/BoostPythonUtil/ListConverter.h"
#include "Python/esys/lsm/geometry/GougeConfigPy.h"
#include "Python/esys/lsm/geometry/GougeConfigPrmsPy.h"
#include "Geometry/VtkXmlWriter.h"

namespace esys
{
  namespace lsm
  {
    GougeConfigPy::GougeConfigPy(const GougeConfigPrmsPy &prms)
      : Inherited(prms)
    {
    }

    void GougeConfigPy::writeVtkXml(const std::string &fileName)
    {
      ParticleDataVisitor particleData;
      visitParticles(particleData);
      visitConnections(particleData);

      VtkXmlWriter vtkWriter;
      vtkWriter.setData(particleData);
      vtkWriter.writeToFile(fileName);
    }

    GougeConfigPy::BBoxVisitor::BBoxVisitor()
      : m_minPt(DBL_MAX), m_maxPt(-DBL_MAX)
    {
    }

    BoundingBoxPy GougeConfigPy::BBoxVisitor::getBBox() const
    {
      return BoundingBoxPy(m_minPt, m_maxPt);
    }

    template <typename TmplParticle>
    void GougeConfigPy::BBoxVisitor::visitSimpleParticle(TmplParticle &particle)
    {
      for (int i = 0; i < 3; i++)
      {
        if ((particle.getPos()[i]-particle.getRad()) < m_minPt[i])
        {
          m_minPt[i] = particle.getPos()[i]-particle.getRad();
        }
        if ((particle.getPos()[i]+particle.getRad()) > m_maxPt[i])
        {
          m_maxPt[i] = particle.getPos()[i]+particle.getRad();
        }
      }
    }

    boost::python::list GougeConfigPy::getConnectionList() const
    {
      return bpu::vectorToList(getConnectionSet());
    }

    boost::python::list GougeConfigPy::getCircDimList() const
    {
      return bpu::vectorToList(getPrms().getPeriodicDimensions());
    }

    BoundingBoxPy GougeConfigPy::getParticleBoundingBox()
    {
      BoundingBoxPy bBox;
      if (getNumParticles() > 0)
      {
        BBoxVisitor bBoxVisitor;
        visitParticles(bBoxVisitor);
        
        bBox = bBoxVisitor.getBBox();
      }
      return bBox;
    }

    BoundingBoxPy GougeConfigPy::getDomainBoundingBox()
    {
      BoundingBoxPy bBox = getParticleBoundingBox();
      Vec3Py minPt       = bBox.getMinPt();
      Vec3Py maxPt       = bBox.getMaxPt();
      for (int i = 0; i < 3; i++)
      {
        if (m_nTablePtr->getPeriodicDimensions()[i])
        {
          minPt[i] = m_nTablePtr->getBBox().getMinPt()[i];
          maxPt[i] = m_nTablePtr->getBBox().getMaxPt()[i];
        }
      }
      if (getPrms().is2d())
      {
        minPt[2] = 0.0;
        maxPt[2] = 0.0;
      }
      return BoundingBoxPy(minPt, maxPt);
    }

    void exportGougeConfig()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<GougeConfigPy, boost::noncopyable>(
        "GougeConfig",
        boost::python::init<const GougeConfigPrmsPy &>()
      )
      .def("generate",                 &GougeConfigPy::generate)
      .def("getNumParticles",          &GougeConfigPy::getNumParticles)
      .def("getNumGrains",             &GougeConfigPy::getNumGrains)
      .def("getParticleCollection",    &GougeConfigPy::getParticleCollection)
      .def("getGrainCollection",       &GougeConfigPy::getGrainCollection)
      .def("getConnectionList",        &GougeConfigPy::getConnectionList)
      .def("write",                    &GougeConfigPy::writeToFile)
      .def("writeVtkXml",              &GougeConfigPy::writeVtkXml)
      .def("getParticleBoundingBox",   &GougeConfigPy::getParticleBoundingBox)
      .def("getParticleBBox",          &GougeConfigPy::getParticleBoundingBox)
      .def("getDomainBoundingBox",     &GougeConfigPy::getDomainBoundingBox)
      .def("getDomainBBox",            &GougeConfigPy::getDomainBoundingBox)
      .def("getCircDimList",           &GougeConfigPy::getCircDimList)
      .def("tagGougeParticles",        &GougeConfigPy::tagGougeParticles)
      .def("tagRndBlockParticles",     &GougeConfigPy::tagRndBlockParticles)
      .def("tagDrivingPlateParticles", &GougeConfigPy::tagDrivingPlateParticles)
      ;
    }
  }
}
