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


#ifndef ESYS_LSMWALLFORCESTOFRICTIONCONVERTER_H
#define ESYS_LSMWALLFORCESTOFRICTIONCONVERTER_H

#include <string>
#include <boost/shared_ptr.hpp>

namespace esys
{
  namespace lsm
  {
    class WallForcesToFrictionConverter
    {
    public:
      WallForcesToFrictionConverter(
        const std::string &wallForcesFile,
        const std::string &instFrictionFile,
        const std::string &avrgFrictionFile,
        int halfWindowSize,
        int wallId1,
        int wallId2,
        int normalDimIndex,
        int shearDimIndex
      );

      void convert();
    private:
      class Impl;
      typedef boost::shared_ptr<Impl> ImplPtr;
      
      ImplPtr m_implPtr;
    };
  }
}

#endif
