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

#ifndef ESYS_LSMWALLPRMSPY_H
#define ESYS_LSMWALLPRMSPY_H

// --- project includes ---
#include "Python/esys/lsm/util/Vec3Py.h"
#include "Model/EWallInteractionGroup.h"
#include "Model/BWallInteractionGroup.h"
#include "Model/SoftBWallInteractionGroup.h"



// --- STL includes ---
#include <string>

using namespace esys::lsm;

namespace esys
{
  namespace lsm
  {
    /*!
      \class NRotElasticWallPrmsPy
      \brief wrapper for CEWallIGP 
      
      $Revision$
      $Date$
    */
    class NRotElasticWallPrmsPy : public  CEWallIGP
    {
    public:
      NRotElasticWallPrmsPy(
        const std::string&,
        const std::string&,
        double
      );
    };

    /*!
      \class NRotBondedWallPrmsPy
      \brief wrapper for CBWallIGP

      $Revision$
      $Date$
    */
    class NRotBondedWallPrmsPy : public CBWallIGP
    { 
    private:
 
    public:
      NRotBondedWallPrmsPy(const std::string&,const std::string&,double,int);
      NRotBondedWallPrmsPy(const std::string&,const std::string&,double,int,int);
    };
    
    /*!
      \class NRotSoftBondedWallPrmsPy
      \brief wrapper for CSoftBWallIGP

      $Revision$
      $Date$
    */
    class NRotSoftBondedWallPrmsPy : public CSoftBWallIGP
    { 
    private:
 
    public:
      NRotSoftBondedWallPrmsPy(const std::string&,const std::string&,double,double,int,int,bool);
      NRotSoftBondedWallPrmsPy(const std::string&,const std::string&,double,double,int,int);
    };

    void exportWallPrms();

  } // namespace lsm
} // namespace esys

#endif // ESYS_LSMWALLPRMSPY_H
