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


#include "Foundation/PathUtil.h"

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/version.hpp>

using namespace boost;

namespace esys
{
  namespace lsm
  {
    void setPathEnv(int argc, char *argv[])
    {
      if (argc > 0)
      {
        setPathEnv(argv[0]);
      }
    }

    /**
     * Function which modifies the PATH environment variable according
     * to the specified executable file. This is a work-around for the
     * SGI MPT mpirun implementation which appears to alter the PATH environment
     * variable of the executed processes.
     */
    void setPathEnv(const std::string &exeName)
    {
      if (exeName.size() > 0)
      {
        std::string origPathValue = "";
        const char *getenvVal = getenv("PATH");
        if (getenvVal != NULL)
        {
          origPathValue = getenvVal;
        }

#if !defined(BOOST_FILESYSTEM_VERSION) || BOOST_FILESYSTEM_VERSION == 2
	filesystem::path exePath = filesystem::system_complete(filesystem::path(exeName, boost::filesystem::native));
	std::string newPathValue = origPathValue + ":" + (exePath.branch_path().native_file_string());
#else 
	filesystem::path exePath = filesystem::system_complete(filesystem::path(exeName));
	std::string newPathValue = origPathValue + ":" + (exePath.branch_path().string());
#endif 
        setenv("PATH", newPathValue.c_str(), 1);
      }
    }
  }
}
