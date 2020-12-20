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


#include "Tools/MacroFrictionCalculator/WallForcesToFrictionConverter.h"
#include "Tools/MacroFrictionCalculator/WallForceReader.h"
#include "Tools/MacroFrictionCalculator/MacroFrictionCalculator.h"
#include "Tools/MacroFrictionCalculator/LinearWindowAverager.h"

#include <fstream>
#include <iomanip>

namespace esys
{
  namespace lsm
  {
    class WallForcesToFrictionConverter::Impl
    {
    public:
      Impl(
        const std::string &wallForcesFile,
        const std::string &instFrictionFile,
        const std::string &avrgFrictionFile,
        int halfWindowSize,
        int wallId1,
        int wallId2,
        int normalDimIndex,
        int shearDimIndex    
      ) :
          m_wallForcesFile(wallForcesFile),
          m_instFrictionFile(instFrictionFile),
          m_avrgFrictionFile(avrgFrictionFile),
          m_halfWindowSize(halfWindowSize),
          m_wallId1(wallId1),
          m_wallId2(wallId2),
          m_normalDimIndex(normalDimIndex),
          m_shearDimIndex(shearDimIndex)
      {
      }

      typedef MacroFrictionCalculator::FrictionVector FrictionVector;
      
      void writeLines(const std::string &fileName, const FrictionVector &vec)
      {
        std::ofstream oStream(fileName.c_str());
        FrictionVector::const_iterator it = vec.begin();
        if (it != vec.end())
        {
          oStream << std::setw(12) << std::setprecision(12) << (*it);
          it++;
        }
        for (; it != vec.end(); it++)
        {
          oStream << "\n" << std::setw(12) << std::setprecision(12) << (*it);
        }
      }

      void convert()
      {
        std::ifstream iStream(m_wallForcesFile.c_str());
        WallForceReader reader(m_wallId1, m_wallId2, iStream);
        
        MacroFrictionCalculator calker(m_normalDimIndex, m_shearDimIndex);
        calker.add(reader);
        
        writeLines(m_instFrictionFile, calker.getFrictionVector());
        
        LinearWindowAverager avgr(calker.getFrictionVector(), m_halfWindowSize, 0, calker.getFrictionVector().size(), 1);
        writeLines(m_avrgFrictionFile, avgr.getAveragedVector());
      }

      const std::string m_wallForcesFile;
      const std::string m_instFrictionFile;
      const std::string m_avrgFrictionFile;
      int               m_halfWindowSize;
      int               m_wallId1;
      int               m_wallId2;
      int               m_normalDimIndex;
      int               m_shearDimIndex;    
    };

    WallForcesToFrictionConverter::WallForcesToFrictionConverter(
      const std::string &wallForcesFile,
      const std::string &instFrictionFile,
      const std::string &avrgFrictionFile,
      int halfWindowSize,
      int wallId1,
      int wallId2,
      int normalDimIndex,
      int shearDimIndex    
    )
      : m_implPtr(
          new Impl(
            wallForcesFile,
            instFrictionFile,
            avrgFrictionFile,
            halfWindowSize,
            wallId1,
            wallId2,
            normalDimIndex,
            shearDimIndex          
          )
        )
    {
    }

    void WallForcesToFrictionConverter::convert()
    {
      m_implPtr->convert();
    }
  }
}
