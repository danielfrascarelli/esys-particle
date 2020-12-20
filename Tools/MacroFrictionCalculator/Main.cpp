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


#include <iostream>
#include <fstream>
#include <stdexcept>

#include "Foundation/StringUtil.h"
#include "Tools/MacroFrictionCalculator/WallForcesToFrictionConverter.h"

using namespace esys::lsm;

int main(int argc, char *argv[])
{
  try
  {
    if (argc > 7) {
      std::string wallForcesFile   = argv[1];
      std::string instFrictionFile = argv[2];
      std::string avrgFrictionFile = argv[3];
      const int halfWindowSize     = StringUtil::to<int>(argv[4]);
      const int wallId1            = StringUtil::to<int>(argv[5]);
      const int wallId2            = StringUtil::to<int>(argv[6]);
      const int normalDimIndex     = StringUtil::to<int>(argv[7]);
      const int shearDimIndex      = StringUtil::to<int>(argv[8]);

      WallForcesToFrictionConverter
        converter(
          wallForcesFile,
          instFrictionFile,
          avrgFrictionFile,
          halfWindowSize,
          wallId1,
          wallId2,
          normalDimIndex,
          shearDimIndex
        );
      converter.convert();
    }
    else {
      std::cerr
        << "Usage: "
        << argv[0] << " wallForcesFile instFrictionFile avrgFrictionFile halfWindowSize wallId1 wallId2 normalDimIndex shearDimIndex" << std::endl
        << "Converts wall-forces record data (id1 F1_x F1_y F1_z id2 F2_x F2_y F2_z ...)" << std::endl
        << "to instantaneuos and effective friction values"
        << std::endl << std::endl
        << "wallForcesFile   - wall-force data read from this file." << std::endl
        << "instFrictionFile - Instantaneous macro friction values written to this file." << std::endl
        << "avrgFrictionFile - Averaged/effective friction values written to this file." << std::endl
        << "halfWindowSize   - Size of the averaging window is 1+(2*halfWindowSize)." << std::endl
        << "wallId1          - integer specifying first wall whose forces are used in instantaneous friction calculation." << std::endl
        << "wallId2          - integer specifying second wall whose forces are used in instantaneous friction calculation." << std::endl
        << "normalDimIndex   - integer {0, 1, 2} specifying the dimension in which the \"normal\"" << std::endl
        << "                   force acts (0 is x dim, 1 is y dim, 2 is z dim)." << std::endl
        << "shearDimIndex    - integer {0, 1, 2} specifying the dimension in which the \"shear\"" << std::endl
        << "                   force acts (0 is x dim, 1 is y dim, 2 is z dim)." << std::endl
        << std::endl << std::endl;
    }
  }
  catch (std::runtime_error &e)
  {
    std::cerr << e.what() << std::endl;
    throw;
  }
  catch (...)
  {
    std::cerr << "Unknown exception." << std::endl;
    throw;
  }
  return 0;
}
