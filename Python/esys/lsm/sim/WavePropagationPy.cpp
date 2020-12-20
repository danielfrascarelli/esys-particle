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

#include <mpi.h>
#include <boost/version.hpp>
#include <boost/python.hpp>
#include <fstream>
#include "Parallel/LatticeMaster.h"
#include "Python/BoostPythonUtil/ListConverter.h"
#include "Python/esys/lsm/sim/WavePropagationPy.h"
#include "Python/esys/lsm/LsmMpiPy.h"

namespace esys
{
  namespace lsm
  {
    class ParticleDataWriter
    {
    public:
      ParticleDataWriter(const std::string &fileName)
        : m_oStream(fileName.c_str())
      {
      }

      template <typename TmplParticle>
      void visitParticle(const TmplParticle &p)
      {
        m_oStream
          << p.getPos()
          << " " << p.getPos()-p.getInitPos()
          << " " << p.getVel()
          << "\n";
      }
      
      template <typename TmplParticle>
      void visitRotParticle(const TmplParticle &p)
      {
        visitParticle(p);
      }

      template <typename TmplParticle>
      void visitRotParticleVi(const TmplParticle &p)
      {
        visitParticle(p);
      }

      template <typename TmplParticle>
      void visitRotThermParticle(const TmplParticle &p)
      {
        visitParticle(p);
      }

    private:
      std::ofstream m_oStream;
    };

    class WavePropagationPy : public LsmMpiPy
    {
    public:
      WavePropagationPy(int numWorkers, const boost::python::list &dimList)
        : LsmMpiPy(numWorkers, dimList)
      {
      }

      void writeParticleDataToFilePyIdList(
          const boost::python::list &idList,
          const std::string &fileName
      )
      {
        ParticleDataWriter visitor(fileName);
        getLatticeMaster().visitParticles(
          bpu::listToVector<int>(idList),
          visitor
        );
      }

      void writeParticleDataToFile(
          const std::string &fileName
      )
      {
        ParticleDataWriter visitor(fileName);
        getLatticeMaster().visitParticles(
          m_idVector,
          visitor
        );
      }

      void setParticleDataIdList(const boost::python::list &idList)
      {
        m_idVector = bpu::listToVector<int>(idList);
      }

      boost::python::list getParticleDataIdList() const
      {
        return bpu::vectorToList(m_idVector);
      }

    private:
      typedef std::vector<int> IdVector;
      IdVector m_idVector;
    };

    void exportWavePropagation()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<WavePropagationPy, boost::python::bases<LsmMpiPy> >(
        "WavePropagation",
        boost::python::init<int, const boost::python::list &>()
      )
      .def(
        "writeParticleDataToFile",
        &WavePropagationPy::writeParticleDataToFile
      )
      .def(
        "writeParticleDataToFile",
        &WavePropagationPy::writeParticleDataToFilePyIdList
      )
      .def(
        "setParticleDataIdList",
        &WavePropagationPy::setParticleDataIdList
      )
      .def(
        "getParticleDataIdList",
        &WavePropagationPy::getParticleDataIdList
      )
      ;
    }
  }
}
