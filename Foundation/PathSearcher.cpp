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


#include "Foundation/PathSearcher.h"
#include "Foundation/StringUtil.h"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/version.hpp>

using namespace boost;

namespace esys
{
  namespace lsm
  {
    PathSearcher::PathSearcher(const std::string &delimitedPathList, const std::string &delim)
      : m_pathVector()
    {
      StringUtil::StringVector searchPaths = StringUtil::splitStrings(delimitedPathList, delim);
      for (StringUtil::StringVector::const_iterator it = searchPaths.begin(); it != searchPaths.end(); it++)
      {
#if !defined(BOOST_FILESYSTEM_VERSION) || BOOST_FILESYSTEM_VERSION == 2
	m_pathVector.push_back(filesystem::path(*it,filesystem::native));
#else 
	m_pathVector.push_back(filesystem::path(*it));
#endif
      }
    }
    
    /*
    bool exists(const std::string &fileName)
    {
      return false;
    }
    */

    filesystem::path PathSearcher::findPath(const std::string &fileName)
    {
      for (PathVector::const_iterator it = m_pathVector.begin(); it != m_pathVector.end(); it++)
      {
        const filesystem::path p = (*it) / fileName;
        if (filesystem::exists(p))
        {
          return p;
        }
      }
      return filesystem::path();
    }

    std::string PathSearcher::find(const std::string &fileName)
    {
      return findPath(fileName).string();
    }
  }
}
