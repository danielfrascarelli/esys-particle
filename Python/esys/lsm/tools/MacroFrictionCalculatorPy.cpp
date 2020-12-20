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
#include <boost/shared_ptr.hpp>
#include "Tools/MacroFrictionCalculator/MacroFrictionCalculator.h"
#include "Tools/MacroFrictionCalculator/LinearWindowAverager.h"
#include "Tools/MacroFrictionCalculator/WallForceReader.h"
#include "Foundation/StlIterator.h"
#include "Foundation/StringUtil.h"
#include "Python/BoostPythonUtil/ListConverter.h"
#include "Python/esys/lsm/tools/MacroFrictionCalculatorPy.h"
#include "Python/esys/lsm/util/Vec3Py.h"

using namespace boost::python;
using namespace esys::lsm;

#include <vector>
#include <fstream>
#include <stdexcept>

namespace esys
{
  namespace lsm
  {
    class WallForcePairPy : public std::pair<Vec3Py, Vec3Py>
    {
    public:
      typedef std::pair<Vec3Py, Vec3Py> WallForcePair;
      WallForcePairPy(const Vec3Py &v1, const Vec3Py &v2)
        : WallForcePair(v1, v2)
      {
      }

      WallForcePairPy(const WallForcePairPy &pair)
        : WallForcePair(pair.first, pair.second)
      {
      }

      WallForcePairPy(const Vec3 &v1, const Vec3 &v2)
        : WallForcePair(v1, v2)
      {
      }

      WallForcePairPy(const WallForceReader::WallForcePair &pair)
        : WallForcePair(pair.first, pair.second)
      {
      }

      WallForcePairPy(const list &l1, const list &l2)
        : WallForcePair()
      {
        const int llen = Vec3Py().len();
        if ((bpu::len(l1) == llen) && (bpu::len(l2) == llen))
        {
          const Vec3Py first1 = 
            Vec3Py(
              extract<double>(l1[0])(),
              extract<double>(l1[1])(),
              extract<double>(l1[2])()
            );
          const Vec3Py second1 = 
            Vec3Py(
              extract<double>(l2[0])(),
              extract<double>(l2[1])(),
              extract<double>(l2[2])()
            );
          first = first1;
          second = second1;
        }
        else
        {
          std::stringstream msg;
          msg 
            << "Lists "
            << extract<std::string>(boost::python::str(l1))()
            << " (len=" << bpu::len(l1) << ")"
            << " and " << extract<std::string>(boost::python::str(l2))()
             << " (len=" << bpu::len(l2) << ")"
            << " do not both have length " << llen;
          throw std::runtime_error(msg.str());
        }
      }

      int len() const
      {
        return 2;
      }

      const Vec3 &getItem(int i) const
      {
        if (i == 0)
        {
          return first;
        }
        else if (i == 1)
        {
          return second;
        }
        throwIndexOutOfRange(i);
        return Vec3::ZERO;
      }

      void setItem(int i, const Vec3 &v)
      {
        if (i == 0)
        {
          first = v;
        }
        else if (i == 1)
        {
          second = v;
        }
        throwIndexOutOfRange(i);
      }

      std::string toString() const
      {
        return 
          extract<std::string>(boost::python::str(first))()
          +
          extract<std::string>(boost::python::str(second))();
      }

    protected:
      void throwIndexOutOfRange(int i) const
      {
        std::stringstream msg;
        msg << "Index i = " << i << " out of range (0,2).";
        throw std::runtime_error(msg.str());        
      }
    };

    class WallForceReaderPy
    {
    public:
      WallForceReaderPy(int wallId1, int wallId2, const std::string &fileName)
        : m_iFStreamPtr(new std::ifstream(fileName.c_str())),
          m_it(wallId1, wallId2, *m_iFStreamPtr)
      {
      }

      bool hasNext() const
      {
        return m_it.hasNext();
      }

      WallForceReader::WallForcePair next()
      {
        return m_it.next();
      }
    private:
      typedef boost::shared_ptr<std::ifstream> IFStreamPtr;
      IFStreamPtr     m_iFStreamPtr;
      WallForceReader m_it;
    };

    class MacroFrictionCalculatorPy
    {
    public:
      MacroFrictionCalculatorPy(
        int normalDimIndex,
        int shearDimIndex
      )
        : m_fricnCalker(normalDimIndex, shearDimIndex)
      {
      }

      class Extractor
      {
      public:
        MacroFrictionCalculator::WallForcePair operator()(object pyOb)
        {
          extract<WallForcePairPy> pair(pyOb);
          if (pair.check())
          {
            return pair();
          }
          else
          {
            std::stringstream msg;
            msg 
              << "Could not extract C++ WallForcePairPy object from list element: "
              << extract<std::string>(str(pyOb))();
            throw runtime_error(msg.str());
          }          
        }

      };

      void addList(const list &wallForceList)
      {
        typedef MacroFrictionCalculator::WallForcePair WallForcePair;
        typedef std::vector<WallForcePair>             WallForceVector;
        
        const WallForceVector 
          wallForceVector(
            bpu::listToVector<MacroFrictionCalculator::WallForcePair, Extractor>(wallForceList)
          );
        m_fricnCalker.add(ForwardConstIterator<WallForceVector>(wallForceVector));
      }

      void addReader(const WallForceReaderPy &wallForceReader)
      {
        m_fricnCalker.add(wallForceReader);
      }

      list getFrictionList() const
      {
        return bpu::vectorToList(m_fricnCalker.getFrictionVector());
      }
      
    private:
      MacroFrictionCalculator m_fricnCalker;
    };

    class LinearWindowAveragerPy
    {
    public:
      LinearWindowAveragerPy(
        const list &valList,
        int halfWindowSize,
        int beginIndex,
        int endIndex,
        int skipSize
      ) : m_linearAverager(
            bpu::listToVector<double>(valList),
            halfWindowSize,
            beginIndex,
            endIndex,
            skipSize
          )
      {
      }

      list getAveragedList()
      {
        return bpu::vectorToList(m_linearAverager.getAveragedVector());
      }
    private:
      LinearWindowAverager m_linearAverager;
    };

    void exportMacroFrictionCalculator()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      class_<WallForceReaderPy>("WallForceReader", init<int,int, const std::string & >());
      class_<WallForcePairPy>("WallForcePair", init<const Vec3Py &, const Vec3Py &>())
        .def(init<const list &,const list &>())
        .def(init<const WallForcePairPy &>())
        .def(self == self)
        .def_readwrite("first",  &WallForcePairPy::first)
        .def_readwrite("second", &WallForcePairPy::second)
        .def("__str__", &WallForcePairPy::toString)
      ;
      class_<MacroFrictionCalculatorPy>("MacroFrictionCalculator", init<int,int>())
        .def("getFrictionList", &MacroFrictionCalculatorPy::getFrictionList)
        .def("add", &MacroFrictionCalculatorPy::addList)
        .def("add", &MacroFrictionCalculatorPy::addReader)
      ;
    
      class_<LinearWindowAveragerPy>("LinearWindowAverager", init<const list &,int,int,int,int>())
        .def("getAveragedList", &LinearWindowAveragerPy::getAveragedList)
      ;
    }
  }
}

